"""Main entrypoint for WDL2CWL."""
import os
from typing import List, Union, Optional, Callable
import WDL  # type: ignore[import]
import cwl_utils.parser.cwl_v1_2 as cwl

from io import StringIO
import textwrap
import re
import sys
import argparse


from ruamel.yaml import scalarstring
from ruamel.yaml.main import YAML


# # WDL-CWL Type Mappings
# wdl_type = {
#     "Array[String]": "string[]",
#     "String": "string",
#     "File": "File",
#     "Int": "int",
#     "Float": "float",
#     "Boolean": "boolean",
# }


# class Converter:
#     @staticmethod
#     def load_wdl_tree(doc: str):
#         wdl_path = os.path.relpath(doc)
#         doc_tree = WDL.load(wdl_path)

#         parser = Converter()

#         if doc_tree.workflow:
#             return parser.load_wdl_objects(doc_tree.workflow)

#         tasks = []
#         for task in doc_tree.tasks:
#             tasks.append(parser.load_wdl_objects(task))

#         return tasks[0]

#     def load_wdl_objects(self, obj: WDL.SourceNode):
#         if isinstance(obj, WDL.Task):
#             return self.load_wdl_task(obj)
#         elif isinstance(obj, WDL.Workflow):
#             return self.load_wdl_workflow(obj)

#     def load_wdl_workflow(self, obj: WDL.Workflow):
#         print(f"Workflow {obj.name} loaded")

#     def load_wdl_task(self, obj: WDL.Task):
#         runtime = obj.runtime

#         cwl_inputs = self.get_cwl_inputs(obj.inputs)
#         base_command = ["bash", "example.sh"]

#         cat_tool = cwl.CommandLineTool(
#             id=obj.name,
#             inputs=cwl_inputs,
#             requirements=None,
#             outputs=[],
#             cwlVersion="v1.2",
#             baseCommand=base_command,
#         )

#         yaml = YAML()
#         yaml.default_flow_style = False
#         yaml.indent = 4
#         yaml.block_seq_indent = 2
#         result_stream = StringIO()
#         cwl_result = cat_tool.save()
#         scalarstring.walk_tree(cwl_result)
#         yaml.dump(cwl_result, result_stream)
#         yaml.dump(cwl_result, sys.stdout)

#         return result_stream.getvalue()

#     def get_cwl_inputs(self, wdl_inputs: List[str]):
#         inputs = []

#         for wdl_input in wdl_inputs:
#             input_name = wdl_input.name
#             input_value = None

#             if isinstance(wdl_input.type, WDL.Type.Array):
#                 input_type = "File"
#                 type_of = [cwl.CommandInputArraySchema(items=input_type, type="array")]
#             elif isinstance(wdl_input.type, WDL.Type.String):
#                 type_of = "string"
#             elif isinstance(wdl_input.type, WDL.Type.Boolean):
#                 type_of = "boolean"
#             elif isinstance(wdl_input.type, WDL.Type.Int):
#                 type_of = "int"
#             else:
#                 type_of = "unknown type"

#             if wdl_input.type.optional:
#                 type_of = [type_of, "null"]

#             if wdl_input.expr is not None:
#                 input_value = wdl_input.expr.literal.value

#             inputs.append(
#                 cwl.CommandInputParameter(
#                     id=input_name, type=type_of, default=input_value
#                 )
#             )

#         return inputs

#     def translate_command(self, expr: WDL.Expr.Base, inputs: List[str]):

#         if expr is None:
#             return None

#         if isinstance(expr, WDL.Expr.Array):
#             return [self.translate_expr(e) for e in expr.items]

#         if isinstance(expr, WDL.Expr.String):
#             return self.translate_command_string(expr)
#         elif isinstance(expr, (WDL.Expr.Int, WDL.Expr.Boolean, WDL.Expr.Float)):
#             return self.literal.value
#         if isinstance(expr, WDL.Expr.Placeholder):
#             return self.translate_expr(expr.expr)

#     def translate_command_string(self, string: WDL.Expr.String):
#         # print("this is the literal", string.literal)
#         # if string.literal is not None:
#         #     return str(string.literal).lstrip("').rstrip('"')

#         # elements = {}
#         # counter = 1
#         # _format = str(string).lstrip('"').rstrip('"')
#         # print("from _format", _format)

#         # for placeholder in string.children:

#         #     print(placeholder)
#         pass
def convert(path: str) -> str:
    """Convert just the bowtie_1 file."""
    wdl_path = path
    doc_tree = WDL.load(wdl_path)
    task = doc_tree.tasks[0]
    inputs = task.inputs
    wdl_first_input = inputs[0]

    cwl_inputs: List[cwl.CommandInputParameter] = []

    # add a function to return the cwl type from wdl_type using isinstance

    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_first_input.name,
            type=[cwl.CommandInputArraySchema(items="File", type="array")],
        )
    )
    wdl_second_input = inputs[1]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_second_input.name,
            type="string",
            default=wdl_second_input.expr.literal.value,
        )
    )
    wdl_third_input = inputs[2]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_third_input.name,
            type=[cwl.CommandInputArraySchema(items="File", type="array")],
        )
    )
    wdl_fourth_input = inputs[3]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_fourth_input.name,
            type="boolean",
            default=wdl_fourth_input.expr.literal.value,
        )
    )
    wdl_fifth_input = inputs[4]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_fifth_input.name,
            type="boolean",
            default=wdl_fifth_input.expr.literal.value,
        )
    )
    wdl_sixth_input = inputs[5]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_sixth_input.name,
            type="boolean",
            default=wdl_sixth_input.expr.literal.value,
        )
    )
    wdl_seventh_input = inputs[6]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_seventh_input.name, type=["int", "null"], default=None
        )
    )
    wdl_eight_input = inputs[7]
    # since wdl_eight_input.type.optional == True then we return an array with null appended
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_eight_input.name, type=["int", "null"], default=None
        )
    )
    wdl_nineth_input = inputs[8]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_nineth_input.name, type=["int", "null"], default=None
        )
    )
    wdl_tenth_input = inputs[9]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_tenth_input.name, type=["string", "null"], default=None
        )
    )
    wdl_eleventh_input = inputs[10]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_eleventh_input.name,
            type="string",
            default=wdl_eleventh_input.expr.literal.value,
        )
    )
    wdl_twelveth_input = inputs[11]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_twelveth_input.name,
            type="int",
            default=wdl_twelveth_input.expr.literal.value,
        )
    )
    wdl_thirteenth_input = inputs[12]
    cwl_inputs.append(
        cwl.CommandInputParameter(
            id=wdl_thirteenth_input.name,
            type="string",
            default=wdl_thirteenth_input.expr.literal.value,
        )
    )

    wdl_outputs = task.outputs
    cwl_outputs = []
    wdl_first_output = wdl_outputs[0]

    cwl_outputs.append(
        cwl.CommandOutputParameter(
            id=wdl_first_output.name,
            type="File",
            outputBinding=cwl.CommandOutputBinding(
                glob="$(inputs.{outputpath})".format(
                    outputpath=wdl_first_output.expr.expr.name
                )
            ),
        )
    )

    requirements: List[cwl.ProcessRequirement] = []

    dockerpull = task.runtime["docker"].expr.referee.expr.literal.value
    requirements.append(cwl.DockerRequirement(dockerPull=dockerpull))
    command = task.command
    cwl_command_str = ""

    wdl_command_part_1 = command.parts[0]
    command_1 = textwrap.dedent(wdl_command_part_1)
    command_1_splitted_p1, command_1_splitted_p2 = command_1.split("$")
    command_1 = command_1_splitted_p1 + "\\$" + command_1_splitted_p2
    cwl_command_str += command_1
    wdl_command_part_2 = command.parts[1]  # This is a wdl.Expr.placeholder object
    command_2_expr = wdl_command_part_2.expr.expr.name
    cwl_command_str += "$(inputs." + command_2_expr + " )"
    wdl_command_part_3 = command.parts[2]
    command_3 = textwrap.dedent(wdl_command_part_3)
    cwl_command_str += command_3
    wdl_command_part_4 = command.parts[3]

    (
        command_4_arg_name,
        command_4_input_value,
    ) = wdl_command_part_4.expr.arguments  # This is a wdl.Expr.apply object
    command_4_arg_name_value = command_4_arg_name.literal.value
    command_4_input_value_name = command_4_input_value.expr.name

    # Since command_4_arg_value.expr.referee.type.optional == True
    cwl_command_str += (
        "$(inputs."
        + command_4_input_value_name
        + ' === null ? "" : '
        + f'"{command_4_arg_name_value}"'
        + " inputs."
        + command_4_input_value_name
        + " )"
    )
    wdl_command_part_5 = command.parts[4]
    command_5 = textwrap.dedent(wdl_command_part_5)
    cwl_command_str += command_5
    wdl_command_part_6 = command.parts[5]

    (
        command_6_arg_name,
        command_6_input_value,
    ) = wdl_command_part_6.expr.arguments  # This is a wdl.Expr.apply object
    command_6_arg_name_value = command_6_arg_name.literal.value
    command_6_input_value_name = command_6_input_value.expr.name

    # Since command_6_arg_value.expr.referee.type.optional == True
    cwl_command_str += (
        "$(inputs."
        + command_6_input_value_name
        + ' === null ? "" : '
        + f'"{command_6_arg_name_value}"'
        + " inputs."
        + command_6_input_value_name
        + " )"
    )

    wdl_command_part_7 = command.parts[6]
    command_7 = textwrap.dedent(wdl_command_part_7)
    cwl_command_str += command_7

    wdl_command_part_8 = command.parts[7]
    (
        command_8_arg_name,
        command_8_input_value,
    ) = wdl_command_part_8.expr.arguments  # This is a wdl.Expr.apply object
    command_8_arg_name_value = command_8_arg_name.literal.value
    command_8_input_value_name = command_8_input_value.expr.name

    # Since command_4_arg_value.expr.referee.type.optional == True
    cwl_command_str += (
        "$(inputs."
        + command_8_input_value_name
        + ' === null ? "" : '
        + f'"{command_8_arg_name_value}"'
        + " inputs."
        + command_8_input_value_name
        + " )"
    )

    wdl_command_part_9 = command.parts[8]
    command_9 = textwrap.dedent(wdl_command_part_9)
    cwl_command_str += command_9
    wdl_command_part_10 = command.parts[9]

    # Since wdl_command_part_10.options = {'true': '--best', 'false': '""'}
    command_10_expr = wdl_command_part_10.expr.expr.name
    cwl_command_str += (
        "$(inputs."
        + command_10_expr
        + " ? "
        + '"{true_value}" : "{false_value}")'.format(
            true_value=wdl_command_part_10.options["true"],
            false_value=wdl_command_part_10.options["false"],
        )
        + ")"
    )
    wdl_command_part_11 = command.parts[10]
    command_11 = textwrap.dedent(wdl_command_part_11)
    cwl_command_str += command_11
    wdl_command_part_12 = command.parts[11]
    command_12_expr = wdl_command_part_12.expr.expr.name
    # Since wdl_command_part_12.options = {'true': '--best', 'false': '""'}
    cwl_command_str += (
        "$(inputs."
        + command_12_expr
        + " ? "
        + '"{true_value}" : "{false_value}")'.format(
            true_value=wdl_command_part_12.options["true"],
            false_value=wdl_command_part_12.options["false"],
        )
        + ")"
    )
    wdl_command_part_13 = command.parts[12]
    command_13 = textwrap.dedent(wdl_command_part_13)
    cwl_command_str += command_13
    wdl_command_part_14 = command.parts[13]
    command_14_expr = wdl_command_part_14.expr.expr.name
    # Since wdl_command_part_14.options = {'true': '--best', 'false': '""
    cwl_command_str += (
        "$(inputs."
        + command_14_expr
        + " ? "
        + '"{true_value}" : "{false_value}")'.format(
            true_value=wdl_command_part_14.options["true"],
            false_value=wdl_command_part_14.options["false"],
        )
        + ")"
    )
    wdl_command_part_15 = command.parts[14]
    command_15 = textwrap.dedent(wdl_command_part_15)
    cwl_command_str += command_15
    wdl_command_part_16 = command.parts[15]

    (
        command_16_arg_name,
        command_16_input_value,
    ) = wdl_command_part_16.expr.arguments  # This is a wdl.Expr.apply object
    command_16_arg_name_value = command_16_arg_name.literal.value
    command_16_input_value_name = command_16_input_value.expr.name

    # Since command_4_arg_value.expr.referee.type.optional == True
    cwl_command_str += (
        "$(inputs."
        + command_16_input_value_name
        + ' === null ? "" : '
        + f'"{command_16_arg_name_value}"'
        + " inputs."
        + command_16_input_value_name
        + " )"
    )

    wdl_command_part_17 = command.parts[16]
    command_17 = textwrap.dedent(wdl_command_part_17)
    cwl_command_str += command_17
    wdl_command_part_18 = command.parts[17]

    (
        command_18_arg_name,
        command_18_input_value,
    ) = wdl_command_part_18.expr.arguments  # This is a wdl.Expr.apply object
    command_18_arg_name_value = command_18_arg_name.literal.value
    command_18_input_value_name = command_18_input_value.expr.name

    # Since command_4_arg_value.expr.referee.type.optional == True
    cwl_command_str += (
        "$(inputs."
        + command_18_input_value_name
        + ' === null ? "" : '
        + f'"{command_18_arg_name_value}"'
        + " inputs."
        + command_18_input_value_name
        + " )"
    )

    # Special cases present...the false_valaue should bome before the true value
    wdl_command_part_19 = command.parts[18]
    command_19_expr = wdl_command_part_19.expr.arguments[0].expr.referee.name
    # Since wdl_command_part_14.options = {'true': '--best', 'false': '""
    cwl_command_str += (
        "$(inputs."
        + command_19_expr
        + " ? "
        + '"{false_value}" : "{true_value}")'.format(
            true_value=wdl_command_part_19.options["true"],
            false_value=wdl_command_part_19.options["false"],
        )
    )

    wdl_command_part_20 = command.parts[19]
    command_20 = textwrap.dedent(wdl_command_part_20)
    cwl_command_str += command_20

    wdl_command_part_21 = command.parts[20]
    wdl_command_part_21_expr = wdl_command_part_21
    wdl_command_part_21_func_name = wdl_command_part_21_expr.expr.function_name
    wdl_command_part_21_arg = wdl_command_part_21_expr.expr.arguments
    (
        wdl_command_part_21_arg_1,
        wdl_command_part_21_arg_2,
        wdl_command_part_21_arg_3,
    ) = wdl_command_part_21_arg

    # since function_name is sub
    (
        wdl_command_part_21_arg_1_input_name,
        wdl_command_part_21_arg_1_index,
    ) = wdl_command_part_21_arg_1.arguments
    command_21 = (
        "$(inputs."
        + wdl_command_part_21_arg_1_input_name.expr.name
        + f"[{wdl_command_part_21_arg_1_index.value}]"
        + f'.replace( "{wdl_command_part_21_arg_2.literal.value}", "{wdl_command_part_21_arg_3.literal.value}") )'
    )
    cwl_command_str += command_21

    wdl_command_part_22 = command.parts[21]
    command_22 = textwrap.dedent(wdl_command_part_22)
    cwl_command_str += command_22

    wdl_command_part_23 = command.parts[22]
    wdl_command_part_23_input_name = wdl_command_part_23.expr.expr.name
    wdl_command_part_23_seperator = wdl_command_part_23.options["sep"]
    command_23 = (
        f"$(inputs.{wdl_command_part_23_input_name}.map("
        + 'function(el) {return el.path}).join("'
        + wdl_command_part_23_seperator
        + '"))'
    )
    cwl_command_str += command_23

    wdl_command_part_24 = command.parts[23]
    command_24 = textwrap.dedent(wdl_command_part_24)
    cwl_command_str += command_24

    wdl_command_part_25 = command.parts[24]
    cwl_command_str += "$(inputs." + wdl_command_part_25.expr.expr.name + ") "

    wdl_command_part_26 = command.parts[25]
    command_26 = textwrap.dedent(wdl_command_part_26)
    cwl_command_str += command_26

    wdl_command_part_27 = command.parts[26]
    cwl_command_str += "$(inputs." + wdl_command_part_27.expr.expr.name + ") "

    wdl_command_part_28 = command.parts[27]
    command_28 = textwrap.dedent(wdl_command_part_28)
    cwl_command_str += command_28

    requirements.append(
        cwl.InitialWorkDirRequirement(
            listing=[cwl.Dirent(entry=cwl_command_str, entryname="example.sh")]
        )
    )
    requirements.append(cwl.InlineJavascriptRequirement())
    requirements.append(cwl.NetworkAccess(networkAccess=True))

    wdl_runtime_cpu = task.runtime["cpu"].expr.name
    ram_min = f"$(inputs.{wdl_runtime_cpu})"
    requirements.append(
        cwl.ResourceRequirement(
            ramMin=ram_min,
        )
    )
    # Resulting cwl output
    base_command = ["bash", "example.sh"]

    cat_tool = cwl.CommandLineTool(
        id=task.name,
        inputs=cwl_inputs,
        requirements=requirements,
        outputs=cwl_outputs,
        cwlVersion="v1.2",
        baseCommand=base_command,
    )

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent = 4
    yaml.block_seq_indent = 2
    result_stream = StringIO()
    cwl_result = cat_tool.save()
    scalarstring.walk_tree(cwl_result)
    yaml.dump(cwl_result, result_stream)
    yaml.dump(cwl_result, sys.stdout)

    return result_stream.getvalue()


def main() -> None:
    """Main entry point."""
    # Command-line parsing.
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow", help="Path to WDL workflow")
    parser.add_argument("-o", "--output", help="Name of resultant CWL file")
    args = parser.parse_args()
    # write to a file in oop_cwl_files
    if args.output:
        with open(args.output, "w") as result:
            result.write(str(convert(args.workflow)))


if __name__ == "__main__":

    main()

    # try:
    #     converted = Converter.load_wdl_tree("wdl2cwl/tests/wdl_files/bowtie_1.wdl")
    #     # converted = Converter.load_wdl_tree(args.workflow)
    # except WDL.Error.SyntaxError as err:
    #     print(err)
    # except WDL.Error.ValidationError as err:
    #     print(err)
    # except WDL.Error.MultipleValidationErrors as err:
    #     for error in err.exceptions:
    #         print(error)
