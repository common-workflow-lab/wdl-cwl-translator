"""Main entrypoint for WDL2CWL."""
import sys
from typing import List, cast

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
}


def get_ram_min(ram_min: str) -> int:
    """
    Get memory requirement.

    Only handles value given in GiB.
    """
    ram_min = ram_min[ram_min.find('"') + 1 : ram_min.find("GiB")]
    return int(float(ram_min.strip()) * 1024)


# need to add sep= and default= for expression placeholders
def get_command(
    command: str,
    unbound: List[str],
    bound: List[str],
    input_types: List[str],
    input_names: List[str],
) -> str:

    index = 0
    new_command = ""
    start_index = 0
    end_index = 0

    # continue till the end of the string
    while index < len(command):

        # if you have ~{
        if command[index] == "~" and command[index + 1] == "{":
            start_index = index + 2

            # while loop to find index of }
            while 1:
                if command[index] == "}":
                    end_index = index
                    break
                else:
                    index += 1

            # sub string containing everything inside ~{ and }
            sub_str = command[start_index:end_index]

            # if sub string has a concatenation
            if "+" in sub_str:
                split_str = sub_str.split("+")

                for i in split_str:
                    if i in input_names:
                        index = input_names.index(i)
                        data_type = input_types[index]  # get the data type of the input

                        if data_type == "File":
                            new_command += "$(inputs." + i + ".path)"
                        else:
                            new_command += "$(inputs." + i + ")"
                    else:
                        new_command += i.replace('"', "")

            # to handle Expression Placeholder Options
            elif ("true" and "false") in sub_str:
                sub = sub_str.split("=")
                true_value = (sub[1].split('"'))[1]
                false_value = (sub[2].split('"'))[1]
                input_name = (sub[2].split('"'))[2]
                append_str = (
                    "$(if ['${inputs."
                    + input_name
                    + "}' eq 'true' ] then echo "
                    + true_value
                    + " else echo "
                    + false_value
                    + " fi )"
                )
                new_command += append_str

            # if sub string has only the input/ variable name
            else:
                data_type = (
                    input_types[input_names.index(sub_str)]
                    if sub_str in input_names
                    else ""
                )
                append_str = ""
                if data_type == "File":
                    append_str = "$(inputs." + sub_str + ".path)"
                else:
                    append_str = "$(inputs." + sub_str + ")"

                new_command = new_command + append_str

            index = end_index + 1
        else:
            new_command += command[index]
            index += 1
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

    input_types: List[str] = []
    input_names: List[str] = []
    input_values: List[str] = []

    for i in ast.task_inputs:
        input_types.append(i[0])
        input_names.append(i[1])
        input_values.append("None")

    for i in ast.task_inputs_bound:
        input_types.append(i[0])
        input_names.append(i[1])
        input_values.append(i[2])

    # returns the entire command including "command{........}"
    raw_command: str = cast(str, ast.task_command)
    raw_command = raw_command[
        raw_command.find("{") + 1 : -1
    ]  # removing the command{} part

    split_command: List[str] = raw_command.strip().split("\n")
    command = ""
    for i in range(0, len(split_command)):
        split_command[i] = split_command[i].strip()
        if "\\" in split_command[i]:
            command += split_command[i]
        elif "#" in split_command[i] or "\n" == split_command[i]:  # skip comments
            continue
        else:
            command += split_command[i] + " "

    command = get_command(
        command, ast.task_inputs, ast.task_inputs_bound, input_types, input_names
    )
    print(command)

    base_command = ["sh", "example.sh"]

    inputs = []
    for i in ast.task_inputs:
        input_type = wdl_type[i[0]]
        input_name = i[1]
        inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))

    requirements: List[cwl.ProcessRequirement] = []

    requirements.append(
        cwl.DockerRequirement(
            dockerPull=ast.task_runtime["docker"].replace('"', ""),
        )
    )

    requirements.append(
        cwl.InitialWorkDirRequirement(
            listing=[cwl.Dirent(entry=command, entryname="example.sh")]
        )
    )

    requirements.append(cwl.InlineJavascriptRequirement())

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
            output_glob = i[2][i[2].find("~{") + 2 : i[2].find("}")]
            output_glob = "$(inputs." + output_glob + ")"
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
        requirements=requirements,
        hints=hints,
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
