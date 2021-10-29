"""Main entrypoint for WDL2CWL."""
import argparse
from argparse import Namespace
import sys
from io import StringIO
from typing import List, cast, Union
from io import StringIO
import textwrap
import re

import regex  # type: ignore

import cwl_utils.parser_v1_2 as cwl
from antlr4 import CommonTokenStream, InputStream  # type: ignore
from ruamel.yaml import scalarstring
from ruamel.yaml.main import YAML

from wdl2cwl.WdlV1_1Lexer import WdlV1_1Lexer
from wdl2cwl.WdlV1_1Parser import WdlV1_1Parser
from wdl2cwl.WdlV1_1ParserVisitor import WdlV1_1ParserVisitor

# WDL-CWL Type Mappings
wdl_type = {
    "Array[String]": "string[]",
    "String": "string",
    "File": "File",
    "Int": "int",
    "Float": "float",
    "Boolean": "boolean",
}

valid_js_identifier = regex.compile(
    r"^(?!(?:do|if|in|for|let|new|try|var|case|else|enum|eval|null|this|true|"
    r"void|with|break|catch|class|const|false|super|throw|while|yield|delete|export|"
    r"import|public|return|static|switch|typeof|default|extends|finally|package|"
    r"private|continue|debugger|function|arguments|interface|protected|implements|"
    r"instanceof)$)(?:[$_\p{ID_Start}])(?:[$_\u200C\u200D\p{ID_Continue}])*$"
)
# ^^ is a combination of https://github.com/tc39/proposal-regexp-unicode-property-escapes#other-examples
# and regex at the bottom of https://stackoverflow.com/a/9392578
# double checked against https://262.ecma-international.org/5.1/#sec-7.6
# eval is not on the official list of reserved words, but it is a built-in function


def inputs(input_name: str) -> str:
    """Produce a consise, valid CWL expr/param reference lookup string for a given input name."""
    if valid_js_identifier.match(input_name):
        return f"inputs.{input_name}"
    return f'inputs["{input_name}"]'


def get_ram_min(ram_min: str) -> int:
    """
    Get memory requirement.

    Only handles value given in GiB.
    """
    unit = " ".join(re.findall("[a-zA-Z]+", ram_min)).replace(" ", "")
    ram_min = ram_min[ram_min.find('"') + 1 : ram_min.find(unit)]
    ram_value = 0
    # Add more units
    if unit == "GiB":
        ram_value = int(float(ram_min.strip()) * 1024)
    return ram_value


def get_ram_min_js(ram_min: str, unit: str) -> str:
    """Get memory requirement for user input."""
    append_str: str = ""
    if unit:
        append_str = '${\nvar unit = "' + unit + '";'
    else:
        append_str = (
            "${\nvar unit = " + inputs(ram_min) + '.match(/[a-zA-Z]+/g).join("");'
        )
    js_str = (
        append_str
        + "\nvar value = parseInt("
        + inputs(ram_min)
        + ".match(/[0-9]+/g));\n"
        + 'var memory = "";\n'
        + 'if(unit==="KiB") memory = value/1024;\n'
        + 'else if(unit==="MiB") memory = value;\n'
        + 'else if(unit==="GiB") memory = value*1024;\n'
        + 'else if(unit==="TiB") memory = value*1024*1024;\n'
        + 'else if(unit==="B") memory = value/(1024*1024);\n'
        + 'else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);\n'
        + 'else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);\n'
        + 'else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);\n'
        + 'else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);\n'
        + "return parseInt(memory);\n}"
    )

    return js_str


def get_command(
    command: str,
    unbound: List[str],
    bound: List[str],
    input_types: List[str],
    input_names: List[str],
    command_check: bool,
) -> str:
    """
    Get command to be used in the bash script.

    Need to add sep= and default= for expression placeholders
    """
    index = 0
    new_command = ""
    start_index = 0
    end_index = 0

    unbound_input_names = []
    for i in unbound:
        unbound_input_names.append(i[1])
    # continue till the end of the string
    while index < len(command):

        # if you have ~{
        if (command[index] == "~" and command[index + 1] == "{") or (
            command[index] == "$" and command[index + 1] == "{" and command_check
        ):
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
                optional_check = False
                optional_inputs = []
                split_str = sub_str.split("+")

                for i in split_str:
                    if '"' not in i:
                        index = input_names.index(i)
                        data_type = input_types[index]

                        if "?" in data_type and i in unbound_input_names:
                            optional_inputs.append(i)
                            optional_check = True

                if optional_check:  # there is one or more optional variable
                    temp_command = ""
                    null_string = ""

                    for i in split_str:
                        if '"' in i:
                            temp_command += i + " + "
                        else:
                            index = input_names.index(i)
                            data_type = input_types[
                                index
                            ]  # get the data type of the input

                            path_str = ""
                            if "File" in data_type:
                                path_str = ".path"

                            if i in optional_inputs:
                                null_string += inputs(i) + " ||"

                            temp_command += inputs(i) + path_str + " + "

                    new_command += (
                        "$("
                        + null_string[:-2]
                        + '=== null ? "" : '
                        + temp_command[:-2]
                        + ")"
                    )
                else:
                    for i in split_str:
                        if '"' in i:
                            new_command += i.replace('"', "")
                        else:
                            index = input_names.index(i)
                            data_type = input_types[index]

                            if data_type != "File":
                                new_command += "$(" + inputs(i) + ")"

            elif "length(" in sub_str:

                if ("true" and "false") in sub_str:
                    true_value = sub_str[
                        sub_str.find("true") + 5 : sub_str.find("false")
                    ].strip()
                    temp = sub_str.split("false=")
                    false_value = temp[1].split(temp[1][0])[1]
                    comparison_expression = temp[1].split(temp[1][0])[2]
                    comparison_expression = comparison_expression[7:]
                    operator = ""
                    for i in comparison_expression:
                        if i in ">=!<":
                            operator += i
                    input_name = comparison_expression.split(operator)[0][:-1]
                    value_to_compare = comparison_expression.split(operator)[1]

                    append_str = (
                        "${if ("
                        + inputs(input_name)
                        + ".length "
                        + operator
                        + " "
                        + value_to_compare
                        + ") {return "
                        + true_value
                        + ";} else {return '"
                        + false_value
                        + "';}}"
                    )
                    new_command += append_str
                else:
                    raise ValueError(
                        "Length function without the if...else keywords is currently not supported"
                    )

            elif ("true" and "false") in sub_str:
                true_value = sub_str[
                    sub_str.find("true") + 5 : sub_str.find("false")
                ].strip()
                temp = sub_str.split("false=")
                false_value = temp[1].split(temp[1][0])[1]
                input_name = temp[1].split(temp[1][0])[2]

                if input_name in input_names:

                    index = input_names.index(input_name)
                    data_type = input_types[index]

                    # true false, when the input is optional and there's no default value assigned
                    if "?" in data_type and input_name in unbound_input_names:

                        append_str = (
                            "$"
                            + inputs(input_name)
                            + ' === null ? "'
                            + false_value
                            + '"'
                            + " : "
                            + true_value
                            + ")"
                        )
                        new_command += append_str

                    else:

                        append_str = (
                            "$("
                            + inputs(input_name)
                            + " ? "
                            + true_value
                            + ' : "'
                            + false_value
                            + '")'
                        )
                        new_command += append_str

                elif "defined(" in input_name:
                    sub_str = input_name[input_name.find("(") + 1 : -1]
                    index = input_names.index(sub_str)
                    data_type = input_types[index]
                    check_str = ""

                    if "Array" in data_type:
                        check_str = ".length === 0 "
                    else:
                        check_str = " === null "

                    append_str = (
                        f"$({inputs(sub_str)}{check_str} ? "
                        f'"{false_value}" : {true_value})'
                    )
                    new_command += append_str

            elif "sep=" in sub_str:
                split_str = sub_str.split(sub_str[4])
                separator = split_str[1]
                input_name = split_str[2]
                index = input_names.index(input_name)
                data_type = input_types[index]

                temp_append_str = ""
                if "Array[File]" in data_type:
                    temp_append_str = ".map(function(el) { return el.path})"

                if input_name in input_names:
                    append_str_sub = (
                        f'({inputs(input_name)}{temp_append_str}.join("{separator}"))'
                    )

                    if "?" in data_type and input_name in unbound_input_names:
                        append_str = f'$({inputs(input_name)} === null ? "" : {(append_str_sub)})'
                    else:
                        append_str = f'$({inputs(input_name)}{temp_append_str}.join("{separator}"))'
                    new_command += append_str

            elif "sub(" in sub_str:
                temp = sub_str.split(",")
                search_index = re.search(r"\[[0-9]\]", temp[0])

                append_str_sub = ""
                if search_index:
                    append_str_sub = (
                        inputs(re.sub(r"\[[0-9]\]", "", temp[0].replace("sub(", "")))
                        + temp[0][temp[0].find("[") :]
                    )
                else:
                    append_str_sub = inputs(temp[0].replace("sub(", ""))

                input_name = ""
                if search_index:
                    input_name = temp[0][: temp[0].find("[")].replace("sub(", "")
                else:
                    input_name = temp[0].replace("sub(", "")

                index = input_names.index(input_name)
                data_type = input_types[index]

                append_str = ""

                if len(temp) == 3:
                    append_str_common = (
                        append_str_sub
                        + ".replace("
                        + temp[1]
                        + ","
                        + temp[2][:-1]
                        + ")"
                    )

                    check_str = ""
                    if "Array" in data_type:
                        check_str = ".length === 0 "

                    if "?" in data_type and input_name in unbound_input_names:
                        append_str = (
                            "$"
                            + inputs(input_name)
                            + check_str
                            + '? "" :'
                            + append_str_common
                            + ")"
                        )
                    else:
                        append_str = "$(" + append_str_common + ")"

                    new_command += append_str
            else:

                data_type = (
                    input_types[input_names.index(sub_str)]
                    if sub_str in input_names
                    else ""
                )

                append_str = ""

                if "?" in data_type and sub_str in unbound_input_names:
                    if "File" in data_type:
                        append_str = f'$({inputs(sub_str)} === null ? "" : {inputs(sub_str)}.path)'
                    else:
                        append_str = (
                            f'$({inputs(sub_str)} === null ? "" : {inputs(sub_str)})'
                        )

                else:

                    if data_type == "File":
                        append_str = f"$({inputs(sub_str)}.path)"
                    else:
                        append_str = "$(inputs." + sub_str + ")"

                new_command = new_command + append_str

            index = end_index + 1
        elif (command[index] == "$" and command[index + 1] == "(") or (
            command[index] == "$" and command[index + 1] == "("
        ):
            new_command += "\\" + command[index]
            index += 1
        else:
            new_command += command[index]
            index += 1
    return new_command


def get_output(expression: str, input_names: List[str]) -> str:
    """Get expression for outputs."""
    output_value = ""

    # for parameter references
    # might have to change, In case there's more than one ~{}
    if "~" in expression and "glob(" not in expression:
        start_index = expression.find("~{")
        end_index = expression.find("}")
        output_value = (
            expression[0:start_index]
            + "$("
            + inputs(expression[start_index + 2 : end_index])
            + ")"
            + expression[end_index + 1 :]
        )
        output_value = output_value.replace('"', "")

    elif "${" in expression:
        start_index = expression.find("${")
        end_index = expression.find("}")
        output_value = (
            expression[0:start_index]
            + "$("
            + inputs(expression[start_index + 2 : end_index])
            + ")"
            + expression[end_index + 1 :]
        )
        output_value = output_value.replace('"', "")
    # For a string concatenation
    elif "+" in expression:
        split_str = expression.split("+")

        for i in split_str:
            if i in input_names:
                output_value += f"$({inputs(i)})"
            else:
                output_value += i

        output_value = output_value.replace('"', "")

    elif "glob(" in expression:
        if "~{" in expression:
            sub_expression = expression.replace("glob(", "")[:-1]
            start_index = sub_expression.find("~{")
            end_index = sub_expression.find("}")
            output_value = (
                sub_expression[0:start_index]
                + "$("
                + inputs(sub_expression[start_index + 2 : end_index])
                + ")"
                + sub_expression[end_index + 1 :]
            )

        else:
            sub_expression = expression.replace("glob(", "")[:-1]
            output_value = sub_expression

        output_value = output_value.replace('"', "")

    elif '"' not in expression:
        if expression in input_names:
            output_value = f"$({inputs(expression)})"
        output_value = output_value

    elif '"' in expression:
        output_value = expression.replace('"', "")

    return output_value


def get_input(
    inputs: List[cwl.CommandInputParameter],
    unbound_input: List[str],
    bound_input: List[str],
) -> List[cwl.CommandInputParameter]:
    """Get bound and unbound inputs."""
    for i in unbound_input:

        input_name = i[1]

        if "Array" in i[0]:
            temp_type = wdl_type[
                i[0][i[0].find("[") + 1 : i[0].find("]")].replace('"', "")
            ]
            input_type = temp_type if "?" not in i[0] else [temp_type, "null"]
            input_name = i[1]

            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name,
                    type=[cwl.CommandInputArraySchema(items=input_type, type="array")],
                )
            )

        else:
            input_type = (
                wdl_type[i[0]]
                if "?" not in i[0]
                else [wdl_type[i[0].replace("?", "")], "null"]
            )

            inputs.append(cwl.CommandInputParameter(id=input_name, type=input_type))

    for i in bound_input:

        input_name = i[1]

        input_type = (
            wdl_type[i[0]]
            if "?" not in i[0]
            else [wdl_type[i[0].replace("?", "")], "null"]
        )

        raw_input_value = i[2].replace('"', "")
        input_value: Union[str, bool, int, float] = ""

        if "boolean" in input_type:
            input_value = bool(raw_input_value.lower() == "true")
        elif "int" in input_type:
            input_value = int(raw_input_value)
        elif "float" in input_type:
            input_value = float(raw_input_value)
        else:
            input_value = raw_input_value

        inputs.append(
            cwl.CommandInputParameter(
                id=input_name,
                type=input_type,
                default=input_value,
            )
        )

    return inputs


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

    # variable to check the syntax of the command { or <<<
    command_check = True
    if "<<<" in raw_command:
        raw_command = raw_command[raw_command.find("<<<") + 3 : -3]
        command_check = False
    else:
        raw_command = raw_command[
            raw_command.find("{") + 1 : -1
        ]  # removing the command{} part

    command = textwrap.dedent(raw_command)

    command = get_command(
        command,
        ast.task_inputs,
        ast.task_inputs_bound,
        input_types,
        input_names,
        command_check,
    )

    base_command = ["bash", "example.sh"]

    cwl_inputs: List[cwl.CommandInputParameter] = []
    cwl_inputs = get_input(cwl_inputs, ast.task_inputs, ast.task_inputs_bound)

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
    requirements.append(cwl.NetworkAccess(networkAccess=True))

    if "memory" in ast.task_runtime:

        ram_min: Union[str, int] = ""

        if '"' in ast.task_runtime["memory"]:
            if "${" in ast.task_runtime["memory"]:
                input_name = ast.task_runtime["memory"][
                    ast.task_runtime["memory"].find("${")
                    + 2 : ast.task_runtime["memory"].find("}")
                ]
                temp = ast.task_runtime["memory"].index("}")  # index of }
                if len(ast.task_runtime["memory"]) != temp + 1:
                    unit = ast.task_runtime["memory"][temp + 1 : -1].strip()
                    if input_name in input_names:
                        ram_min = get_ram_min_js(input_name, unit)
            else:
                ram_min = get_ram_min(ast.task_runtime["memory"])
        else:
            ram_min = get_ram_min_js(ast.task_runtime["memory"], "")

        requirements.append(
            cwl.ResourceRequirement(
                ramMin=ram_min,
            )
        )

    if "time_minutes" in ast.task_runtime:

        time_minutes: Union[str, int] = ""
        if '"' not in ast.task_runtime["time_minutes"]:
            if ast.task_runtime["time_minutes"] in input_names:
                time_minutes = f"$({inputs(ast.task_runtime['time_minutes'])} * 60)"
            elif ast.task_runtime["time_minutes"].isnumeric():
                time_minutes = int(ast.task_runtime["time_minutes"]) * 60
        else:
            print(
                f"---- WARNING: SKIPPING REQUIREMENT time_minutes in  {ast.task_runtime} ----",
                file=sys.stderr,
            )

        requirements.append(
            cwl.ToolTimeLimit(
                timelimit=time_minutes,
            )
        )

    if "cpu" in ast.task_runtime:

        cpu: Union[str, int] = ""
        if '"' not in ast.task_runtime["cpu"]:
            if ast.task_runtime["cpu"] in input_names:
                cpu = f"$({inputs(ast.task_runtime['cpu'])})"
            elif ast.task_runtime["cpu"].isnumeric():
                cpu = int(ast.task_runtime["cpu"])
            elif "+" in ast.task_runtime["cpu"]:
                temp = ast.task_runtime["cpu"].split("+")
                append_str = "$("
                for i in range(0, len(temp)):
                    if temp[i] in input_names:
                        append_str += inputs(temp[i])
                    elif temp[i].isnumeric():
                        append_str += temp[i]
                    if i != len(temp) - 1:
                        append_str += " + "
                append_str += ")"
                cpu = append_str

        requirements.append(
            cwl.ResourceRequirement(
                coresMin=cpu,
            )
        )

    outputs = []

    for i in ast.task_outputs:
        output_name = i[1]
        output_glob = get_output(i[2], input_names)

        if "Array" in i[0]:
            # output_type = ""
            temp_type = wdl_type[
                i[0][i[0].find("[") + 1 : i[0].find("]")].replace('"', "")
            ]
            output_type = temp_type if "?" not in i[0] else [temp_type, "null"]
            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type=[
                        cwl.CommandOutputArraySchema(items=output_type, type="array")
                    ],
                    outputBinding=cwl.CommandOutputBinding(glob=output_glob),
                )
            )
        elif "read_string(" in i[2]:
            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type="string",
                    outputBinding=cwl.CommandOutputBinding(
                        glob=output_glob.replace("read_string(", "")[:-1],
                        loadContents=True,
                        outputEval=r"$(self[0].contents.replace(/[\r\n]+$/, ''))",
                    ),
                )
            )
        elif i[2] == "stdout()":
            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type="stdout",
                )
            )
        else:
            output_type = (
                wdl_type[i[0]]
                if "?" not in i[0]
                else [wdl_type[i[0].replace("?", "")], "null"]
            )

            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type=output_type,
                    outputBinding=cwl.CommandOutputBinding(glob=output_glob),
                )
            )

    cat_tool = cwl.CommandLineTool(
        id=ast.task_name,
        inputs=cwl_inputs,
        requirements=requirements if requirements else None,
        outputs=outputs if outputs else [],
        cwlVersion="v1.2",
        baseCommand=base_command,
    )

    # implemented runtime requirements
    runtime_requirements = ["docker", "memory", "time_minutes", "cpu"]

    for i in ast.task_runtime:
        if i not in runtime_requirements:
            print("----WARNING: SKIPPING REQUIREMENT " + i + "----", file=sys.stderr)

    if ast.task_parameter_meta_check:
        print("----WARNING: SKIPPING PARAMETER_META----", file=sys.stderr)

    if ast.task_meta_check:
        print("----WARNING: SKIPPING META----", file=sys.stderr)

    if len(ast.task_variables) > 0:
        for a in ast.task_variables:
            print(
                "----WARNING: SKIPPING VARIABLE " + str(a[1]) + "----", file=sys.stderr
            )

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent = 4
    yaml.block_seq_indent = 2
    result_stream = StringIO()
    cwl_result = cat_tool.save()
    scalarstring.walk_tree(cwl_result)
    # ^ converts multiline strings to nice multiline YAML
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
