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
    ram_min = ram_min[ram_min.find('"') + 1 : ram_min.find("GiB")]
    return int(float(ram_min.strip()) * 1024)

#Does not handle cases like ~{true="-2 " false="" twoPassMode}  yet
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
            new_command=new_command+append_str
            index=(end_index+1)
        else:
            new_command+=command[index]
            index+=1
    return new_command

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
        else:
            command+=raw_command[i]+" "

    command = get_command(command,ast.task_inputs,ast.task_inputs_bound)

    base_command = ["sh", "example.sh"]

    inputs = []
    for i in ast.task_inputs:
        input_type = wdl_type[i[0].replace("?", "")]
        input_name = i[1]
        inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))

    for i in ast.task_inputs_bound:
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
    if ast.task_runtime:
        requirements.append(
            cwl.DockerRequirement(
                dockerPull=ast.task_runtime["docker"].replace('"', "")
            )
        )

    requirements.append(
        cwl.InitialWorkDirRequirement(listing=[cwl.Dirent(entry=command,entryname="example.sh")]))

    hints = []
    if ast.task_runtime:
        hints.append(
            cwl.ResourceRequirement(
                ramMin=get_ram_min(ast.task_runtime["memory"]),
            )
        )

    outputs = []

    for i in ast.task_outputs:
        output_type = wdl_type[i[0]]
        output_name = i[1]
        output_glob = ""
        if "~" in i[2]:
            start_index = i[2].find("~{")
            end_index = i[2].find("}")
            output_glob = (
                i[2][0:start_index]
                + "$(inputs."
                + i[2][start_index + 2 : end_index]
                + ")"
                + i[2][end_index + 1 :]
            )
            output_glob = output_glob.replace('"', "")

        else:
            output_glob = i[2]

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

    if "preemptible" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT PREEMPTIBLE----")

    if "disks" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT DISKS----")

    if len(ast.task_variables) > 0:
        for a in ast.task_variables:
            print("----WARNING: SKIPPING VARIABLE " + str(a[1]) + "----")

    return cast(str, yaml.main.round_trip_dump(cat_tool.save()))


if __name__ == "__main__":
    main(sys.argv[1:])
