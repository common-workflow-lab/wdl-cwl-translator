"""Main entrypoint for WDL2CWL."""
import argparse
from argparse import Namespace
import sys
from io import StringIO
from typing import List, cast, Any
from io import StringIO
import textwrap
import re

import cwl_utils.parser_v1_2 as cwl
from antlr4 import CommonTokenStream, InputStream  # type: ignore
from ruamel.yaml import scalarstring
from ruamel.yaml.main import YAML

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
    unit = " ".join(re.findall("[a-zA-Z]+", ram_min))
    ram_min = ram_min[ram_min.find('"') + 1 : ram_min.find(unit)]
    return int(float(ram_min.strip()) * 1024)


def get_command(
    command: str,
    unbound: List[str],
    bound: List[str],
    input_types: List[str],
    input_names: List[str],
) -> str:
    """
    Get command to be used in the bash script.

    Need to add sep= and default= for expression placeholders
    """
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

            # if sub string has only the input/ variable name
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


def convert(workflow: str) -> str:
    """Generate a CWL object to match "cat-tool.cwl"."""
    f = open(workflow)
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

    command = textwrap.dedent(raw_command)

    command = get_command(
        command, ast.task_inputs, ast.task_inputs_bound, input_types, input_names
    )

    base_command = ["sh", "example.sh"]

    inputs = []
    for i in ast.task_inputs:
        input_type = wdl_type[i[0]]
        input_name = i[1]
        inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))

    requirements: List[cwl.ProcessRequirement] = []

    if "docker" in ast.task_runtime:
        dockerPull = ""

        if '"' not in ast.task_runtime["docker"]:
            # only searching for value in bound inputs.
            # value could be a bound declaration which is not handled
            for sublist in ast.task_inputs_bound:
                if ast.task_runtime["docker"] in sublist[1]:
                    dockerPull = sublist[2]

        else:
            dockerPull = ast.task_runtime["docker"]

        requirements.append(
            cwl.DockerRequirement(dockerPull=dockerPull.replace('"', ""))
        )

    requirements.append(
        cwl.InitialWorkDirRequirement(
            listing=[cwl.Dirent(entry=command, entryname="example.sh")]
        )
    )

    requirements.append(cwl.InlineJavascriptRequirement())

    hints: List[cwl.ProcessRequirement] = []

    if "memory" in ast.task_runtime:
        memory: Any = ""

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

    if "time_minutes" in ast.task_runtime:

        time_minutes = ""
        if '"' not in ast.task_runtime["time_minutes"]:
            time_minutes = "$(inputs."+ast.task_runtime["time_minutes"]+"* 60)"

        requirements.append(
            cwl.ToolTimeLimit(
                timelimit=time_minutes.replace('"',""),
            )
        )

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
        requirements=requirements if requirements else None,
        hints=hints if hints else None,
        outputs=outputs if outputs else None,
        cwlVersion="v1.2",
        baseCommand=base_command,
    )

    if ast.task_parameter_meta_check:
        print("----WARNING: SKIPPING PARAMETER_META----")

    if "preemptible" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT PREEMPTIBLE----")

    if "disks" in ast.task_runtime:
        print("----WARNING: SKIPPING REQUIREMENT DISKS----")

    if len(ast.task_variables) > 0:
        for a in ast.task_variables:
            print("----WARNING: SKIPPING VARIABLE " + str(a[1]) + "----")

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent = 4
    yaml.block_seq_indent = 2
    result_stream = StringIO()
    cwl_result = cat_tool.save()
    scalarstring.walk_tree(cwl_result)
    # ^ converts multine line strings to nice multiline YAML
    yaml.dump(cwl_result, result_stream)
    yaml.dump(cwl_result, sys.stdout)

    return result_stream.getvalue()


def main() -> None:
    """Command-line parsing."""
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow", help="Path to WDL workflow")
    parser.add_argument("-o", "--output", help="Name of resultant CWL file")
    args = parser.parse_args()

    if args.output:
        with open(args.output, "w") as result:
            result.write(str(convert(args.workflow)))


if __name__ == "__main__":
    main()
