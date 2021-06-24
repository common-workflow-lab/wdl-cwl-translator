"""Main entrypoint for WDL2CWL."""
import sys
from typing import List, cast
import re

import cwl_utils.parser_v1_2 as cwl
from antlr4 import CommonTokenStream, InputStream  # type: ignore
from ruamel import yaml

from wdl2cwl.WdlV1_1Lexer import WdlV1_1Lexer
from wdl2cwl.WdlV1_1Parser import WdlV1_1Parser
from wdl2cwl.WdlV1_1ParserVisitor import WdlV1_1ParserVisitor

# WDL-CWL Type Mappings
wdl_type = {
    "String": "string",
    "File": "File",
    "Int": "int",
    "Float": "float",
    "Boolean": "boolean",
    "?": "null",
}


def get_ram_min(ram_min: str) -> int:
    """
    Get memory requirement.

    Only handles value given in GiB.
    """
    unit = " ".join(re.findall("[a-zA-Z]+", ram_min))
    ram_min = ram_min[ram_min.find('"') + 1 : ram_min.find(unit)]
    return int(float(ram_min.strip()) * 1024)

def get_command(command,unbound,bound,input_types,input_names):
    
    index = 0
    new_command = ""
    start_index = 0
    end_index = 0
   
   #continue till the end of the string
    while index<len(command):
        
        #if you have ~{
        if command[index] is "~" and command[index+1] is "{":
            start_index = index+2

            #while loop to find index of }
            while 1:
                if command[index] is "}":
                    end_index = index
                    break
                else:
                    index+=1

            #sub string containing everything inside ~{ and }
            sub_str = command[start_index:end_index]
            
            #if sub string has a concatenation
            if "+" in sub_str:
                split_str = sub_str.split("+")
                
                for i in split_str:
                    if i in input_names:
                        index = input_names.index(i)
                        data_type = input_types[index] #get the data type of the input

                        if data_type == "File":
                            new_command+="$(inputs."+i+".path)"
                        else:
                            new_command+="$(inputs."+i+")"
                    else:
                        new_command+=i.replace('"','')
           
           #to handle Expression Placeholder Options
           #might have to change later
            elif ("true" and "false") in sub_str:
                sub= sub_str.split("=")
                true_value =  (sub[1].split('"'))[1]
                false_value = (sub[2].split('"'))[1]
                input_name = (sub[2].split('"'))[2]
                append_str = "$(if ['${inputs."+input_name+"}' eq 'true' ] then echo "+ true_value+" else echo "+ false_value+" fi )"
                new_command+=append_str
                    
            #if sub string has only the input/ variable name
            else:               
                data_type = input_types[input_names.index(sub_str)] if sub_str in input_names else ""
                append_str = ""
                if data_type == "File":
                    append_str = "$(inputs."+sub_str+".path)"
                else:
                    append_str = "$(inputs."+sub_str+")"

                new_command=new_command+append_str

            index=(end_index+1)
        else:
            new_command+=command[index]
            index+=1
    return new_command

def get_workflow(ast):

    inputs = []
    outputs = []
    steps = []
    cat_tool = cwl.Workflow(
        inputs = inputs,
        outputs = outputs,
        steps = steps,
        cwlVersion= "v1.2",
    )
    
    return cast(str, yaml.main.round_trip_dump(cat_tool.save()))

def get_expression_placeholder(expression, input_names):
    expression_value = ""
    #for parameter references
    if "~" in expression:
        start_index = expression.find("~{")
        end_index = expression.find("}")
        expression_value = (
            expression[0:start_index]
             + "$(inputs."
            + expression[start_index + 2 : end_index]
            + ")"
            + expression[end_index + 1 :]
        )
        expression_value = expression_value.replace('"', "")

    #For a string concatenation
    elif "+" in expression:
        split_str = expression.split("+")

        for i in split_str:
            if i in input_names:
                expression_value+="$(inputs."+i+")"
            else:
                expression_value+=i

        expression_value = expression_value.replace('"', "")
        
    else:
        expression_value = expression

    return expression_value

def main(argv: List[str]) -> str:
    """Generate a CWL object to match "cat-tool.cwl"."""
    f = open(argv[0])
    text = InputStream(f.read())
    lexer = WdlV1_1Lexer(text)
    stream = CommonTokenStream(lexer)
    parser = WdlV1_1Parser(stream)
    tree = parser.document()  # type: ignore

    # create WdlV1_1ParserVisitor object and return all inputs, outputs, etc
    ast = WdlV1_1ParserVisitor()  # type: ignore
    ast.walk_tree(tree)  # type: ignore

    input_types = []
    input_names = []
    input_values = []

    for i in ast.task_inputs:
        input_types.append(i[0])
        input_names.append(i[1])
        input_values.append(None)

    for i in ast.task_inputs_bound:
        input_types.append(i[0])
        input_names.append(i[1])
        input_values.append(i[2])

    # returns the entire command including "command{........}"
    raw_command: str = cast(str, ast.task_command)
    raw_command = raw_command[
        raw_command.find("{") + 1 : -1
    ]  # removing the command{} part

    raw_command = raw_command.strip().split('\n')
    command = ""
    for i in range(0,len(raw_command)):
        raw_command[i] = raw_command[i].strip()
        if "\\" in raw_command[i]:
            command+=raw_command[i]
        elif "#" in raw_command[i] or "\n" is raw_command[i]: #skip comments
            continue
        else:
            command+=raw_command[i]+" "
        
    command = get_command(command,ast.task_inputs,ast.task_inputs_bound,input_types,input_names)
    print(command)

    base_command = ["sh", "example.sh"]

    inputs = []
    for i in ast.task_inputs:           
        if "Array" in i[0]:
            input_type = i[0][i[0].find('[')+1:-1].replace('"',"")
            input_name = i[1]
            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name, 
                    type=[cwl.CommandInputArraySchema(
                        items=input_type, type="array"
                    )]
                )
            )
        else:
            input_type = wdl_type[i[0].replace("?", "")]
            input_name = i[1]
            inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))


    for i in ast.task_inputs_bound:
        if "Array" in i[0]:
            input_type = i[0][i[0].find('[')+1:-1].replace('"',"")
            input_name = i[1]
            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name, 
                    type=[cwl.CommandInputArraySchema
                        (items=input_type, type="array")], 
                    default=i[2]
                )
            )
        else:
            input_type = wdl_type[i[0].replace("?", "")]
            input_name = i[1]
            input_expression = i[2].replace('"', "")
            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name,
                    type=input_type,
                    default=input_expression,
                )
            )

    requirements = []
    if "docker" in ast.task_runtime:
        dockerPull = ""

        if '"' not in ast.task_runtime["docker"]:
            for sublist in ast.task_inputs_bound:
                if ast.task_runtime["docker"] in sublist[1]:
                    dockerPull = sublist[2]

        elif '"' in ast.task_runtime["docker"] and "~{" in ast.task_runtime["docker"]:
            start_index = ast.task_runtime["docker"].find("~{")
            end_index = ast.task_runtime["docker"].find("}")
            sub_str = ast.task_runtime["docker"][start_index+2:end_index]

            if sub_str in input_names:
                index = input_names.index(sub_str)
                dockerPull = input_values[index]

        else:
            dockerPull = ast.task_runtime["docker"]
        
        requirements.append(
            cwl.DockerRequirement(
                dockerPull=dockerPull.replace('"', "") 
            )
        )

    requirements.append(
        cwl.InitialWorkDirRequirement(listing=[cwl.Dirent(entry=command,entryname="example.sh")]))

    requirements.append(
        cwl.InlineJavascriptRequirement()
    )

    hints = []
    if "memory" in ast.task_runtime:
        memory = ""

        if '"' not in ast.task_runtime["memory"]:
            for sublist in ast.task_inputs_bound:
                if sublist[1] in ast.task_runtime["memory"]:
                    memory = sublist[2]
        else:
            memory = ast.task_runtime["memory"]

        hints.append(
            cwl.ResourceRequirement(
                ramMin=get_ram_min(memory),
            )
        )

    outputs = []

    for i in ast.task_outputs:

        if "Array" in i[0]:
            output_name = i[1]
            output_type = i[0][i[0].find('[')+1:-1].replace('"',"")
            output_glob = get_expression_placeholder(i[2], input_names)
            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type=[cwl.CommandOutputArraySchema(
                        items=output_type, type="array"
                    )],
                    outputBinding=cwl.CommandOutputBinding(glob=output_glob),
                )
            )
        else:
            output_type = wdl_type[i[0]]
            output_name = i[1]
            output_glob = get_expression_placeholder(i[2],input_names)

            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type=output_type,
                    outputBinding=cwl.CommandOutputBinding(glob=output_glob),
                )
            )

    cat_tool = cwl.CommandLineTool(
        id=ast.task_name,
        inputs=inputs,
        requirements=requirements if requirements else None,
        hints=hints if hints else None,
        outputs=outputs,
        cwlVersion="v1.2",
        baseCommand=base_command,
    )

    if ast.task_parameter_meta_check:
        print("----WARNING: SKIPPING PARAMETER_META----")

    if ast.task_meta_check:
        print("----WARNING: SKIPPING META----")
        
    if "preemptible" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT PREEMPTIBLE----")

    if "disks" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT DISKS----")

    if "time_minutes" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT TIME_MINUTES----")

    if len(ast.task_variables) > 0:
        for a in ast.task_variables:
            print("----WARNING: SKIPPING VARIABLE " + str(a[1]) + "----")

    with open('result.cwl', 'w') as result:
        result.write(yaml.main.round_trip_dump(cat_tool.save()))

    return cast(str, yaml.main.round_trip_dump(cat_tool.save()))


if __name__ == "__main__":
    main(sys.argv[1:])
