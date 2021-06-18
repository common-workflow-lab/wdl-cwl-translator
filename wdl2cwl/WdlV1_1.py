import sys
from antlr4 import *
from wdl2cwl.WdlV1_1Lexer import WdlV1_1Lexer
from wdl2cwl.WdlV1_1Parser import WdlV1_1Parser
from wdl2cwl.WdlV1_1ParserVisitor import WdlV1_1ParserVisitor
import cwl_utils.parser_v1_2 as cwl

from ruamel import yaml

###### WDL-CWL Type Mappings #########
wdl_type = {
        "String": "string",
        "File": "File",
        "Int": "int",
        "Float": "float",
        "Boolean": "boolean",
        "?": "null"
    }

#get memory requirement
#only to handle values given in GiB
def get_ram_min(ram_min):
    ram_min = ram_min[ram_min.find("\"")+1:ram_min.find("GiB")]
    return int(float(ram_min.strip())*1024)

def get_command(command,unbound,bound):
    index = 0
    new_command = ""
    start_index = 0
    end_index = 0

    input_types = []
    input_names = []
    
    for i in unbound:
        input_types.append(i[0])
        input_names.append(i[1])

    for i in bound:
        input_types.append(i[0])
        input_names.append(i[1])
   
    while index<len(command):
        if command[index]!="\n":
            new_command+=command[index]
        if command[index] is "~" and command[index+1] is "{":
            start_index = index+2
            while 1:
                if command[index] is "}":
                    end_index = index
                    index+=1
                    break
                else:
                    index+=1
            sub_str = command[start_index:end_index]
            data_type = input_types[input_names.index(sub_str)] if sub_str in input_names else ""
            append_str = ""
            if data_type == "File":
                append_str = "$(inputs."+sub_str+".path)"
            else:
                append_str = "$(inputs."+sub_str+")"
            new_command=new_command[:-1]+append_str
        index+=1
    return new_command
    


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
    command = get_command(command,ast.task_inputs,ast.task_inputs_bound)
    command = command.strip("\n")

    '''command = command.strip().split("\\") #split by '\'

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
    '''

    base_command = ["sh", "example.sh"]

    inputs = []
    for i in ast.task_inputs:
        input_type = wdl_type[i[0]] if '?' not in i[0] else [ wdl_type[i[0].replace('?','')], wdl_type["?"] ]
        input_name = i[1]
        inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))  

    for i in ast.task_inputs_bound:
        input_type = wdl_type[i[0]] if '?' not in i[0] else [ wdl_type[i[0].replace('?','')], wdl_type["?"].replace('"','') ]
        input_name = i[1]
        input_expression = i[2].replace('"','')
        inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type, 
        default=input_expression,))  

    requirements = []
    if ast.task_runtime:
        requirements.append(cwl.DockerRequirement(
                dockerPull=ast.task_runtime["docker"].replace('"', '')
            ))
    
    requirements.append(
        cwl.InitialWorkDirRequirement(listing=cwl.Dirent(entry=command,entryname="example.sh")))

    hints = []
    if ast.task_runtime:
        hints.append(cwl.ResourceRequirement(
                ramMin=get_ram_min(ast.task_runtime["memory"]),
            ))   

    outputs = []

    for i in ast.task_outputs:
        output_type = wdl_type[i[0]]
        output_name = i[1]
        output_glob = ""
        if "~" in i[2]:
            start_index = i[2].find("~{")
            end_index = i[2].find("}")
            output_glob = i[2][0:start_index]+"$(inputs."+i[2][start_index+2:end_index]+")"+i[2][end_index+1:]
            output_glob = output_glob.replace('"','')
            
        else:
            output_glob = i[2]

        outputs.append(cwl.CommandOutputParameter(
            id=output_name,
            type=output_type,
            outputBinding=cwl.CommandOutputBinding(glob=output_glob),
        ))

    '''arguments = []

    for i in command_arguments:
        arguments.append(cwl.CommandLineBinding(valueFrom=i))'''
    

    cat_tool = cwl.CommandLineTool(
        id=ast.task_name,
        inputs=inputs,
        requirements=requirements if requirements else None,
        hints=hints if hints else None,
        outputs=outputs,
        cwlVersion="v1.0",
        baseCommand=base_command,
        #arguments=arguments,
        
    )

    if "preemptible" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT PREEMPTIBLE----")

    if "disks" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT DISKS----")

    if(len(ast.task_variables)>0):
        for a in ast.task_variables:
            print("----WARNING: SKIPPING VARIABLE "+str(a[1])+"----")

    with open('result.cwl', 'w') as result:
        result.write(yaml.main.round_trip_dump(cat_tool.save()))
        
    return yaml.main.round_trip_dump(cat_tool.save())
    
if __name__ == "__main__":
    main(sys.argv[1:])