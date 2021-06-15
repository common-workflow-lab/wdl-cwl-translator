import sys
from antlr4 import *
from WdlV1_1Lexer import WdlV1_1Lexer
from WdlV1_1Parser import WdlV1_1Parser
from WdlV1_1ParserVisitor import WdlV1_1ParserVisitor
import cwl_utils.parser_v1_2 as cwl

from ruamel import yaml

###### WDL-CWL Type Mappings #########
wdl_type = {
        "String": "string",
        "File": "File",
        "Int": "int",
        "Float": "float",
        "Boolean": "boolean",
    }

#get memory requirement
#only to handle values given in GiB
def get_ram_min(ram_min):
    ram_min = ram_min[ram_min.find("\"")+1:ram_min.find("GiB")]
    return int(float(ram_min.strip())*1024)

def main(argv) -> None:
    """Generate a CWL object to match "cat-tool.cwl"."""
    
    f = open(argv[0], "r")
    text = InputStream(f.read())
    lexer = WdlV1_1Lexer(text)
    stream = CommonTokenStream(lexer)
    parser = WdlV1_1Parser(stream)
    tree = parser.document()

    #create WdlV1_1ParserVisitor object and return all inputs, outputs, etc
    ast = WdlV1_1ParserVisitor()
    ast.walk_tree(tree)

    #returns the entire command including "command{........}"
    command = ast.task_command
    command = command[command.find("{")+1:-1] #removing the command{} part
    command = command.strip().split("\\") #split by '\'

    base_command = ""
    command_arguments = []

    for a in command:
        #if it contains = , it's taken as an argument
        if "=" in a or "~" in a:
            command_arguments.append(a.strip())
        #else add to base command
        else:
            base_command+=a

    #split the base command by spaces
    base_command = base_command.split()

    ## get command arguments
    index = 0
    for i in command_arguments:
        if "~" in i:
            parameter_reference = i[i.find("~{")+2:i.find("}")] #finding the name of the parameter
            sub_str = i.strip().split("~")
            if "INPUT" not in sub_str[0]:
                command_arguments[index] = sub_str[0]+"$(inputs."+parameter_reference+")"
            else:
                command_arguments[index] = sub_str[0]+"$(inputs."+parameter_reference+".path"")"
        index+=1

    inputs = []
    for i in ast.task_inputs:
        input_type = wdl_type[i[0]]
        input_name = i[1]
        inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))  

    docker_requirement = [
        cwl.DockerRequirement(
            dockerPull=ast.task_runtime["docker"].replace('"', ''),
        ),
    ]

    hints = [
        cwl.ResourceRequirement(
            ramMin=get_ram_min(ast.task_runtime["memory"]),
        )
    ]

    outputs = []

    for i in ast.task_outputs:
        output_type = wdl_type[i[0]]
        output_name = i[1]
        output_glob = ""
        if "~" in i[2]:
            output_glob = i[2][i[2].find("~{")+2:i[2].find("}")]
            output_glob = "$(inputs."+output_glob+")"
        else:
            output_glob = i[2]

        outputs.append(cwl.CommandOutputParameter(
            id=output_name,
            type=output_type,
            outputBinding=cwl.CommandOutputBinding(glob=output_glob),
        ))

    arguments = []

    for i in command_arguments:
        arguments.append(cwl.CommandLineBinding(valueFrom=i))
    

    cat_tool = cwl.CommandLineTool(
        id=ast.task_name,
        inputs=inputs,
        requirements=docker_requirement,
        hints=hints,
        outputs=outputs,
        cwlVersion="v1.0",
        baseCommand=base_command,
        arguments=arguments,
        
    )

    if "preemptible" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT PREEMPTIBLE----")

    if "disks" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT DISKS----")

    if(len(ast.task_variables)>0):
        for a in ast.task_variables:
            print("----WARNING: SKIPPING VARIABLE "+str(a[1])+"----")

    return yaml.main.round_trip_dump(cat_tool.save())
    
if __name__ == "__main__":
    main(sys.argv[1:])
