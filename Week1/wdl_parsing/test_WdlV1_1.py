f = open("Qc.wdl", "r")

import sys
from antlr4 import *
from WdlV1_1Lexer import WdlV1_1Lexer
from WdlV1_1Parser import WdlV1_1Parser
from WdlV1_1ParserVisitor import WdlV1_1ParserVisitor
import cwl_utils.parser_v1_2 as cwl

text = InputStream(f.read())
lexer = WdlV1_1Lexer(text)
stream = CommonTokenStream(lexer)
parser = WdlV1_1Parser(stream)
tree = parser.document()

#create WdlV1_1ParserVisitor object and return all inputs, outputs, etc
ast = WdlV1_1ParserVisitor()
ast = ast.walk_tree(tree)

#get input type and name
def get_inputs(ast):
    print('TASK INPUTS')
    for a in ast.task_inputs:
        print(a[0], " ", a[1])

#get runtime requirements
def get_runtime(ast):
    print('TASK RUNTIME')
    for a in ast.task_runtime:
        print(a[0]," ",a[1])

#get output type, name, expression
def get_outputs(ast):
    print('TASK OUTPUT')
    for a in ast.task_outputs:
        print(a[0]," ",a[1]," ",a[2])

#get command
def get_command(ast):
    print('TASK COMMAND')
    print(ast.task_command)

#get_inputs(ast)
#get_runtime(ast)
#get_outputs(ast)
#get_command(ast)

#returns the entire command including "command{........}"
command = ast.task_command
command = command[command.find("{")+1:-1] #removing the command{} part
command = command.strip().split("\\") #split by '\'

base_command = ""
command_arguments = []

for a in command:
    #if it contains = , it's taken as an argument
    if "=" in a:
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
        command_arguments[index] = sub_str[0]+"$(inputs."+parameter_reference+")"
    index+=1

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

import cwl_utils.parser_v1_2 as cwl

from ruamel import yaml

def main() -> None:
    """Generate a CWL object to match "cat-tool.cwl"."""
    
    inputs = []
    for i in ast.task_inputs:
        input_type = wdl_type[i[0]]
        input_name = i[1]
        inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))  

    docker_requirement = [
        cwl.DockerRequirement(
            dockerPull=ast.task_runtime["docker"],
        )
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
        id="CollectQualityYieldMetrics",
        inputs=inputs,
        requirements=docker_requirement,
        hints=hints,
        outputs=outputs,
        cwlVersion="v1.0",
        baseCommand=base_command,
        arguments=arguments,
    )

    print(yaml.main.round_trip_dump(cat_tool.save()))

if __name__ == "__main__":
    main()


'''python create_cwl_from_objects.py > result.cwl
cwltool --validate result.cwl'''

