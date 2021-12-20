"""Main entrypoint for WDL2CWL."""
import os
from typing import List, Union, Optional, Callable, cast
import WDL
import cwl_utils.parser.cwl_v1_2 as cwl

from io import StringIO
import textwrap
import sys
import argparse


from ruamel.yaml import scalarstring
from ruamel.yaml.main import YAML


class Converter:
    """Object that handles WDL Workflows and task conversion to CWL."""

    @staticmethod
    def load_wdl_tree(doc: str) -> str:
        """Load WDL file, instantiate Converter class and loads the WDL document tree."""
        wdl_path = os.path.relpath(doc)
        doc_tree = WDL.load(wdl_path)

        parser = Converter()

        if doc_tree.workflow:
            return parser.load_wdl_objects(doc_tree.workflow)

        tasks = []
        for task in doc_tree.tasks:
            tasks.append(parser.load_wdl_objects(task))

        return tasks[0]

    def load_wdl_objects(self, obj: Union[WDL.Tree.Task, WDL.Tree.Workflow]) -> str:
        """Load a WDL SourceNode obj and returns either a Task or a Workflow."""
        if isinstance(obj, WDL.Tree.Task):
            return self.load_wdl_task(obj)
        raise Exception(f"Unimplemented type: {type(obj)}: {obj}")

    #     elif isinstance(obj, WDL.Workflow):
    #         return self.load_wdl_workflow(obj)

    # def load_wdl_workflow(self, obj: WDL.Workflow):
    #     print(f"Workflow {obj.name} loaded")
    #     pass

    def load_wdl_task(self, obj: WDL.Tree.Task) -> str:
        """Load task and convert to CWL."""
        cwl_inputs = self.get_cwl_inputs(obj.inputs)
        cwl_outputs = self.get_cwl_outputs(obj.outputs)
        runtime_docker = obj.runtime["docker"]
        if not isinstance(runtime_docker, WDL.Expr.Placeholder):
            raise Exception(
                f"Unsupport docker runtime type: {type(runtime_docker)}: {runtime_docker}"
            )
        docker_requirement = self.get_cwl_docker_requirements(runtime_docker)
        cwl_command_str = self.get_cwl_command_requirements(obj.command.parts)
        base_command = ["bash", "example.sh"]
        requirements: List[cwl.ProcessRequirement] = [
            docker_requirement,
            cwl_command_str,
        ]
        requirements.append(cwl.InlineJavascriptRequirement())
        requirements.append(cwl.NetworkAccess(networkAccess=True))
        cpu_requirement = self.get_cpu_requirement(obj.runtime["cpu"])
        requirements.append(cpu_requirement)

        cat_tool = cwl.CommandLineTool(
            id=obj.name,
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

    def get_cpu_requirement(
        self, cpu_runtime: WDL.Expr.Base
    ) -> cwl.ResourceRequirement:
        """Translate WDL Runtime CPU requirement to CWL Resource Requirement."""
        if not isinstance(cpu_runtime, WDL.Expr.Placeholder):
            raise Exception(f"Unhandled type: {type(cpu_runtime)}: {cpu_runtime}")
        cpu_runtime_name = cast(WDL.Expr.Ident, cpu_runtime.expr).name
        ram_min = f"$(inputs.{cpu_runtime_name})"
        return cwl.ResourceRequirement(ramMin=ram_min)

    def get_cwl_docker_requirements(
        self, wdl_docker: WDL.Expr.Placeholder
    ) -> cwl.ProcessRequirement:
        """Translate WDL Runtime Docker requirements to CWL Docker Requirement."""
        dockerpull_expr = wdl_docker.expr
        if dockerpull_expr is None or not isinstance(dockerpull_expr, WDL.Expr.Ident):
            raise Exception(
                f"Unsupported type: {type(dockerpull_expr)}: {dockerpull_expr}"
            )
        dockerpull_referee = dockerpull_expr.referee
        if dockerpull_referee is None:
            raise Exception(f"Unsupported type: {type(dockerpull_referee)}")
        dockerpull = dockerpull_referee.expr.literal.value
        return cwl.DockerRequirement(dockerPull=dockerpull)

    def get_cwl_command_requirements(
        self, wdl_commands: List[Union[str, WDL.Expr.Placeholder]]
    ) -> cwl.InitialWorkDirRequirement:
        """Translate WDL commands into CWL Initial WorkDir REquirement."""
        command_str: str = ""
        for wdl_command in wdl_commands:
            if isinstance(wdl_command, str):
                command_str += self.translate_wdl_str(wdl_command)
            elif isinstance(wdl_command, WDL.Expr.Placeholder):
                command_str += self.translate_wdl_placeholder(wdl_command)
        return cwl.InitialWorkDirRequirement(
            listing=[cwl.Dirent(entry=command_str, entryname="example.sh")]
        )

    def translate_wdl_placeholder(self, wdl_placeholder: WDL.Expr.Placeholder) -> str:
        """Translate WDL Expr Placeholder to a valid CWL command string."""
        cwl_command_str = ""

        options = wdl_placeholder.options
        if options:
            if "true" in options:
                true_value = options["true"]
                false_value = options["false"]
            elif "sep" in options:
                seperator = options["sep"]
        if isinstance(wdl_placeholder.expr, WDL.Expr.Get):
            nested_expr = wdl_placeholder.expr.expr
            if nested_expr is None or not isinstance(nested_expr, WDL.Expr.Ident):
                raise Exception(f"Unsupported type: {type(nested_expr)}")
            placeholder_name = nested_expr.name
            if not options:
                cwl_command_str = "$(inputs." + placeholder_name + ")"
            elif "true" in options:
                cwl_command_str = (
                    "$(inputs."
                    + placeholder_name
                    + " ? "
                    + f'"{true_value}" : "{false_value}")'
                )
            elif "sep" in options:
                cwl_command_str = (
                    f"$(inputs.{placeholder_name}.map("
                    + 'function(el) {return el.path}).join("'
                    + seperator
                    + '"))'
                )
        elif isinstance(wdl_placeholder.expr, WDL.Expr.Apply):
            function_name = wdl_placeholder.expr.function_name
            expr_arguments = wdl_placeholder.expr.arguments

            if function_name == "defined":
                arg = expr_arguments[0]
                if not isinstance(arg, WDL.Expr.Placeholder):
                    raise Exception(f"Unsupported type: {type(arg)}: {arg}")
                arg_expr = arg.expr
                if not isinstance(arg_expr, WDL.Expr.Ident):
                    raise Exception(f"Unsupported type: {type(arg_expr)}: {arg_expr}")
                arg_referee = arg_expr.referee
                if not isinstance(arg_expr, WDL.Tree.Decl):
                    raise Exception(f"Unsupported type: {type(arg_expr)}: {arg_expr}")
                arg_referee_name = arg_referee.name
                cwl_command_str = (
                    "$(inputs."
                    + arg_referee_name
                    + " ? "
                    + f'"{false_value}" : "{true_value}")'
                )
            elif function_name == "_interpolation_add":
                arg_name_raw, arg_value_raw = expr_arguments
                arg_name_literal = arg_name_raw.literal
                if arg_name_literal is None or not hasattr(arg_name_literal, "value"):
                    raise Exception(f"Unsupported type: {type(arg_name_literal)}")
                if not isinstance(arg_value_raw, WDL.Expr.Placeholder):
                    raise Exception(f"Unsupported type: {type(arg_value_raw)}")
                arg_value_expr = arg_value_raw.expr
                if not isinstance(arg_value_expr, WDL.Expr.Ident):
                    raise Exception(f"Unsupported type: {type(arg_value_expr)}")
                arg_name, arg_value = arg_name_literal.value, arg_value_expr.name
                if wdl_placeholder.expr.type.optional:
                    cwl_command_str = (
                        "$(inputs."
                        + arg_value
                        + ' === null ? "" : '
                        + f'"{arg_name}"'
                        + " inputs."
                        + arg_value
                        + " )"
                    )
                else:
                    cwl_command_str = f"{arg_name} $(inputs.{arg_value})"
            elif function_name == "sub":
                wdl_apply_object_arg, arg_string, arg_substitute = expr_arguments
                apply_input, index_to_sub = wdl_apply_object_arg.arguments
                apply_input_name, index_to_sub_value = (
                    apply_input.expr.name,
                    index_to_sub.value,
                )
                cwl_command_str = (
                    "$(inputs."
                    + apply_input_name
                    + f"[{index_to_sub_value}]"
                    + f'.replace("{arg_string.literal.value}", "{arg_substitute.literal.value}") )'
                )
        return cwl_command_str

    def translate_wdl_str(self, wdl_command: str) -> str:
        """Translate WDL string command to CWL Process requirement string."""
        first_newline = wdl_command.find("\n")
        command_str = wdl_command[:first_newline] + textwrap.dedent(
            wdl_command[first_newline:]
        )

        if "$" in command_str:
            splitted_1, splitted_2 = command_str.split("$")
            command_str = splitted_1 + "\\$" + splitted_2

        return command_str

    def get_cwl_inputs(
        self, wdl_inputs: Optional[List[WDL.Tree.Decl]]
    ) -> List[cwl.CommandInputParameter]:
        """Convert WDL inputs into CWL inputs and return a list of CWL Command Input Paramenters."""
        inputs: List[cwl.CommandInputParameter] = []

        for wdl_input in wdl_inputs:
            input_name = wdl_input.name
            input_value = None
            type_of: Union[str, cwl.CommandInputArraySchema]

            if isinstance(wdl_input.type, WDL.Type.Array):
                input_type = "File"
                type_of = cwl.CommandInputArraySchema(items=input_type, type="array")

            elif isinstance(wdl_input.type, WDL.Type.String):
                type_of = "string"
            elif isinstance(wdl_input.type, WDL.Type.Boolean):
                type_of = "boolean"
            elif isinstance(wdl_input.type, WDL.Type.Int):
                type_of = "int"
            else:
                type_of = "unknown type"

            if wdl_input.type.optional:
                final_type_of: Union[
                    List[Union[str, cwl.CommandInputArraySchema]],
                    str,
                    cwl.CommandInputArraySchema,
                ] = [type_of, "null"]
            else:
                final_type_of = type_of

            if wdl_input.expr is not None:
                input_value = wdl_input.expr.literal.value

            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name, type=final_type_of, default=input_value
                )
            )

        return inputs

    def get_cwl_outputs(
        self, wdl_outputs: List[WDL.Tree.Decl]
    ) -> List[cwl.CommandOutputParameter]:
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


def main() -> None:
    """Entry point."""
    # Command-line parsing.
    # parser = argparse.ArgumentParser()
    # parser.add_argument("workflow", help="Path to WDL workflow")
    # parser.add_argument("-o", "--output", help="Name of resultant CWL file")
    # args = parser.parse_args()

    # # write to a file in oop_cwl_files
    # if args.output:
    #     with open(args.output, "w") as result:
    # result.write(str(Converter.load_wdl_tree("wdl2cwl/tests/wdl_files/bcftools_stats.wdl"))) #missing args.workflow)

    Converter.load_wdl_tree("wdl2cwl/tests/wdl_files/bcftools_stats.wdl")


if __name__ == "__main__":

    main()
