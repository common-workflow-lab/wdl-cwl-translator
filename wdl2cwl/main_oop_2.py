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


class Converter:
    """An object that handles WDL Workflows and task conversion to CWL."""
    @staticmethod
    def load_wdl_tree(doc: str):
        """A static method that loads the WDL file and loads the WDL document tree."""
        wdl_path = os.path.relpath(doc)
        doc_tree = WDL.load(wdl_path)

        parser = Converter()

        if doc_tree.workflow:
            return parser.load_wdl_objects(doc_tree.workflow)

        tasks = []
        for task in doc_tree.tasks:
            tasks.append(parser.load_wdl_objects(task))

        return tasks[0]

    def load_wdl_objects(self, obj: WDL.SourceNode):
        """Loads a WDL SourceNode obj and returns either a Task or a Workflow."""
        if isinstance(obj, WDL.Task):
            return self.load_wdl_task(obj)
        elif isinstance(obj, WDL.Workflow):
            return self.load_wdl_workflow(obj)

    # def load_wdl_workflow(self, obj: WDL.Workflow):
    #     print(f"Workflow {obj.name} loaded")
    #     pass

    def load_wdl_task(self, obj: WDL.Task):
        """Load task and convert to CWL."""

        cwl_inputs = self.get_cwl_inputs(obj.inputs)
        cwl_outputs = self.get_cwl_outputs(obj.outputs)
        docker_requirement = self.get_cwl_docker_requirements(obj.runtime["docker"])
        test_part_11 = obj.command.parts[11]
        test_part_18 = obj.command.parts[18]
        test_part_20 = obj.command.parts[20]
        test_part_22 = obj.command.parts[22]
        test_placeholder_with_args = obj.command.parts[3]
        cwl_command_str = self.get_cwl_command_requirements(obj.command.parts)
        base_command = ["bash", "example.sh"]

        cat_tool = cwl.CommandLineTool(
            id=obj.name,
            inputs=cwl_inputs,
            requirements=docker_requirement + cwl_command_str,
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

    def get_cwl_docker_requirements(self, wdl_docker: WDL.Expr):
        requirements: List[cwl.ProcessRequirement] = []
        dockerpull = wdl_docker.expr.referee.expr.literal.value
        requirements.append(cwl.DockerRequirement(dockerPull=dockerpull))
        return requirements

    def get_cwl_command_requirements(self, wdl_commands: List[str]) -> List[cwl.InitialWorkDirRequirement]:
        """Translate WDL commands into CWL Initial WorkDir REquirement."""
        command_str: str = ""
        for wdl_command in wdl_commands:
            if isinstance(wdl_command, str):
                command_str += self.translate_wdl_str(wdl_command)
            elif isinstance(wdl_command, WDL.Expr.Placeholder):
                command_str += self.translate_wdl_placeholder(wdl_command)
        return [
            cwl.InitialWorkDirRequirement(
                listing=[cwl.Dirent(entry=command_str, entryname="example.sh")]
            )
        ]
    def translate_wdl_placeholder(self, wdl_placeholder) -> str:
        cwl_command_str = ""
        if isinstance(wdl_placeholder.expr, WDL.Expr.Get):
            placeholder_name = wdl_placeholder.expr.expr.name
            if wdl_placeholder.options:
                if "true" in wdl_placeholder.options:
                    cwl_command_str = (
            "$(inputs."
            + placeholder_name
            + " ? "
            + '"{true_value}" : "{false_value}")'.format(
                true_value=wdl_placeholder.options["true"],
                false_value=wdl_placeholder.options["false"],
            )
            + ")"
        )
            else:
                cwl_command_str = "$(inputs." + wdl_placeholder.expr.expr.name + ")"
        elif isinstance(wdl_placeholder.expr, WDL.Expr.Apply):
            if len(wdl_placeholder.expr.arguments) == 1:
                arg_name = wdl_placeholder.expr.arguments[0].expr.referee.name
                cwl_command_str = (
                    "$(inputs."
                    + arg_name
                    + " ? "
                    + '"{false_value}" : "{true_value}")'.format(
                        true_value=wdl_placeholder.options["true"],
                        false_value=wdl_placeholder.options["false"],
                    )
                )
            if len(wdl_placeholder.expr.arguments) == 2:
                arg_name, arg_value = wdl_placeholder.expr.arguments
                arg_name, arg_value = arg_name.literal.value, arg_value.expr.name
                if wdl_placeholder.expr.type.optional:
                    cwl_command_str = (
                        "$(inputs."
                        + arg_value
                        + ' === null ? "" : '
                        + f'"{arg_name}"'
                        + " inputs."
                        + arg_name
                        + " )"
                    )
        return cwl_command_str

    def translate_wdl_str(self, wdl_command: str) -> str:
        """Translate WDL string command to CWL Process requirement string."""
        wdl_command.replace("\\n", "\n")
        command_str = textwrap.dedent(wdl_command)
        if "$" in command_str:
            splitted_1, splitted_2 = command_str.split("$")
            command_str = splitted_1 + "\\$" + splitted_2


        return command_str

    def get_cwl_inputs(self, wdl_inputs: List[str]) -> List[cwl.CommandInputParameter]:
        """Convert WDL inputs into CWL inputs and return a list of CWL Command Input Paramenters."""
        inputs: List[cwl.CommandInputParameter] = []

        for wdl_input in wdl_inputs:
            input_name = wdl_input.name
            input_value = None

            if isinstance(wdl_input.type, WDL.Type.Array):
                input_type = "File"
                type_of = [cwl.CommandInputArraySchema(items=input_type, type="array")]
            elif isinstance(wdl_input.type, WDL.Type.String):
                type_of = "string"
            elif isinstance(wdl_input.type, WDL.Type.Boolean):
                type_of = "boolean"
            elif isinstance(wdl_input.type, WDL.Type.Int):
                type_of = "int"
            else:
                type_of = "unknown type"

            if wdl_input.type.optional:
                type_of = [type_of, "null"]

            if wdl_input.expr is not None:
                input_value = wdl_input.expr.literal.value

            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name, type=type_of, default=input_value
                )
            )

        return inputs

    def get_cwl_outputs(self, wdl_outputs: List[str]) -> List[cwl.CommandOutputParameter]:
        """Convert WDL outputs into CWL outputs and return a list of CWL Command Output Parameters."""
        outputs: List[cwl.CommandOutputParameter] = []

        for wdl_output in wdl_outputs:
            output_name = wdl_output.name
            if isinstance(wdl_output.type, WDL.Type.File):
                type_of = "File"

            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type=type_of,
                    outputBinding=cwl.CommandOutputBinding(
                        glob="$(inputs.{outputpath})".format(
                            outputpath=wdl_output.expr.expr.name
                        )
                    ),
                )
            )
        return outputs

    # def translate_command(self, expr: WDL.Expr.Base, inputs: List[str]):

    #     if expr is None:
    #         return None

    #     # if isinstance(expr, WDL.Expr.Array):
    #     #     return [self.translate_expr(e) for e in expr.items]

    #     if isinstance(expr, WDL.Expr.String):
    #         return self.translate_command_string(expr)
    #     elif isinstance(expr, (WDL.Expr.Int, WDL.Expr.Boolean, WDL.Expr.Float)):
    #         return self.literal.value
    #     # if isinstance(expr, WDL.Expr.Placeholder):
    #     #     return self.translate_expr(expr.expr)

    # def translate_command_string(self, string: WDL.Expr.String):

    #     pass


def main() -> None:
    """Entry point."""
    # # Command-line parsing.
    # parser = argparse.ArgumentParser()
    # parser.add_argument("workflow", help="Path to WDL workflow")
    # parser.add_argument("-o", "--output", help="Name of resultant CWL file")
    # args = parser.parse_args()
    # # write to a file in oop_cwl_files
    # if args.output:
    # with open(args.output, "w") as result:
    #     result.write(str(
    Converter.load_wdl_tree("wdl2cwl/tests/wdl_files/bowtie_1.wdl")


if __name__ == "__main__":

    main()
