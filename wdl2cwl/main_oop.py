"""Main entrypoint for WDL2CWL."""
import os
from typing import List, Union, Optional, Callable, cast, Any
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
        if not isinstance(runtime_docker, WDL.Expr.Get):
            raise Exception(
                f"Unsupported docker runtime type: {type(runtime_docker)}: {runtime_docker}"
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
        if "memory" in obj.runtime and isinstance(obj.runtime["memory"], WDL.Expr.Get):
            memory_requirement = self.get_memory_requirement(obj.runtime["memory"])
        else:
            memory_requirement = None
        requirements.append(
            cwl.ResourceRequirement(
                coresMin=cpu_requirement,
                ramMin=memory_requirement,
            )
        )

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

    def get_memory_requirement(self, memory_runtime: WDL.Expr.Get) -> str:
        """Translate WDL Runtime Memory requirement to CWL Resource Requirement."""
        ram_min = ""
        if isinstance(memory_runtime.expr, WDL.Expr.Ident):
            expr_referee = memory_runtime.expr.referee
            if expr_referee:
                expr_referee_name = str(expr_referee.name)
                ram_min = self.get_ram_min_js(expr_referee_name, "")
        return ram_min

    def get_ram_min_js(self, ram_min_ref_name: str, unit: str) -> str:
        """Get memory requirement for user input."""
        append_str: str = ""
        if unit:
            append_str = '${\nvar unit = "' + unit + '";'
        else:
            append_str = (
                "${\nvar unit = inputs."
                + ram_min_ref_name
                + '.match(/[a-zA-Z]+/g).join("");'
            )
        js_str = (
            append_str
            + "\nvar value = parseInt(inputs."
            + ram_min_ref_name
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

    def get_cpu_requirement(self, cpu_runtime: WDL.Expr.Base) -> str:
        """Translate WDL Runtime CPU requirement to CWL Resource Requirement."""
        if isinstance(cpu_runtime, WDL.Expr.Get):
            cpu_runtime_name = cast(WDL.Expr.Ident, cpu_runtime.expr).name
            cores_min = f"$(inputs.{cpu_runtime_name})"
        elif isinstance(cpu_runtime, WDL.Expr.Apply):
            ref_function = cpu_runtime.function_name
            ref_arguments = cpu_runtime.arguments
            if ref_function == "_add":
                first_arg, second_arg = ref_arguments
                second_arg_value = self.get_wdl_literal(second_arg.literal)
                first_arg_expr_name = cast(
                    WDL.Expr.Ident, cast(WDL.Expr.Get, first_arg).expr
                ).name
                cores_min = f"$(inputs.{first_arg_expr_name} + {second_arg_value})"

        else:
            raise Exception(f"Unhandled type: {type(cpu_runtime)}: {cpu_runtime}")
        return cores_min

    def get_cwl_docker_requirements(
        self, wdl_docker: WDL.Expr.Get
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
                if nested_expr.referee and isinstance(
                    nested_expr.referee.expr, WDL.Expr.Apply
                ):
                    reference_expr = nested_expr.referee
                    ref_expr_to_apply = reference_expr.expr
                    ref_function = ref_expr_to_apply.function_name
                    ref_arguments = ref_expr_to_apply.arguments
                    if ref_function == "_add":
                        # return true value for a javascript tenary
                        first_arg, second_arg = ref_arguments
                        second_arg_value = self.get_wdl_literal(second_arg.literal)
                        first_arg_fun_name = first_arg.function_name
                        if first_arg_fun_name == "basename":
                            only_argument = first_arg.arguments[0]
                            only_argument_expr_name = only_argument.expr.name
                        true_tenary = f"inputs.{only_argument_expr_name}.{first_arg_fun_name} + '{second_arg_value}'"

                    cwl_command_str = (
                        " $(inputs."
                        + reference_expr.name
                        + " === null ?"
                        + f" {true_tenary} : inputs.{reference_expr.name})"
                    )
                else:
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
                if not isinstance(arg, WDL.Expr.Get):
                    raise Exception(f"Unsupported type: {type(arg)}: {arg}")
                arg_expr = arg.expr
                if not isinstance(arg_expr, WDL.Expr.Ident):
                    raise Exception(f"Unsupported type: {type(arg_expr)}: {arg_expr}")
                arg_referee = arg_expr.referee
                if not isinstance(arg_referee, WDL.Tree.Decl):
                    raise Exception(
                        f"Unsupported type: {type(arg_referee)}: {arg_referee}"
                    )
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
                if not isinstance(arg_value_raw, WDL.Expr.Get):
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
                        + ")"
                    )
                else:
                    cwl_command_str = f"{arg_name} $(inputs.{arg_value})"
            elif function_name == "sub":
                wdl_apply_object_arg, arg_string, arg_substitute = expr_arguments
                if not isinstance(wdl_apply_object_arg, WDL.Expr.Apply):
                    raise Exception(f"Unsupported type: {type(wdl_apply_object_arg)}")
                apply_input, index_to_sub = wdl_apply_object_arg.arguments
                if not isinstance(apply_input, WDL.Expr.Get):
                    raise Exception(f"Unsupported type: {type(apply_input)}")
                apply_input_name = apply_input.expr.name  # type: ignore[attr-defined]
                if not hasattr(index_to_sub, "value"):
                    raise Exception(f"{type(index_to_sub)} has no attribute: 'value")
                index_to_sub_value = index_to_sub.value  # type: ignore[attr-defined]
                cwl_command_str = (
                    "$(inputs."
                    + apply_input_name
                    + f"[{index_to_sub_value}]"
                    + f'.replace("{self.get_wdl_literal(arg_string.literal)}", "{self.get_wdl_literal(arg_substitute.literal)}"))'
                )
        return cwl_command_str

    def get_wdl_literal(self, wdl_expr: Optional[WDL.Value.Base]) -> Any:
        """Extract Literal value from WDL expr."""
        if wdl_expr is None or not hasattr(wdl_expr, "value"):
            raise Exception(f"{type(wdl_expr)} has not attribute 'value'")
        return wdl_expr.value

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

        if not wdl_inputs:
            return inputs

        for wdl_input in wdl_inputs:
            input_name = wdl_input.name
            input_value = None
            type_of: Union[str, cwl.CommandInputArraySchema]

            if isinstance(wdl_input.type, WDL.Type.Array):
                input_type = ""
                array_items_type = wdl_input.type.item_type
                if isinstance(array_items_type, WDL.Type.File):
                    input_type = "File"
                elif isinstance(array_items_type, WDL.Type.String):
                    input_type = "string"
                elif isinstance(array_items_type, WDL.Type.Boolean):
                    input_type = "boolean"
                elif isinstance(array_items_type, WDL.Type.Int):
                    input_type = "int"
                else:
                    raise Exception(
                        f"Array of item_type = {type(array_items_type)}: not supported"
                    )
                type_of = cwl.CommandInputArraySchema(items=input_type, type="array")

            elif isinstance(wdl_input.type, WDL.Type.File):
                type_of = "File"
            elif isinstance(wdl_input.type, WDL.Type.String):
                type_of = "string"
            elif isinstance(wdl_input.type, WDL.Type.Boolean):
                type_of = "boolean"
            elif isinstance(wdl_input.type, WDL.Type.Int):
                type_of = "int"
            else:
                type_of = "unknown type"

            if wdl_input.type.optional or isinstance(wdl_input.expr, WDL.Expr.Apply):
                final_type_of: Union[
                    List[Union[str, cwl.CommandInputArraySchema]],
                    str,
                    cwl.CommandInputArraySchema,
                ] = [type_of, "null"]
            else:
                final_type_of = type_of

            if wdl_input.expr is not None:
                if isinstance(wdl_input.expr, WDL.Expr.Apply):
                    input_value = None
                else:
                    literal = wdl_input.expr.literal
                    if not literal or not hasattr(literal, "value"):
                        raise Exception(f'{type(literal)} has no attribute "value"')
                    input_value = literal.value

            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name, type=final_type_of, default=input_value
                )
            )

        return inputs

    def get_expr_name(self, wdl_expr: WDL.Tree.Decl) -> str:
        """Extract name from WDL expr."""
        if wdl_expr is None or not hasattr(wdl_expr, "name"):
            raise Exception(f"{type(wdl_expr)} has not attribute 'name'")
        return wdl_expr.name

    def get_cwl_outputs(
        self, wdl_outputs: List[WDL.Tree.Decl]
    ) -> List[cwl.CommandOutputParameter]:
        """Convert WDL outputs into CWL outputs and return a list of CWL Command Output Parameters."""
        outputs: List[cwl.CommandOutputParameter] = []

        if not wdl_outputs:
            return outputs

        for wdl_output in wdl_outputs:
            output_name = wdl_output.name
            if isinstance(wdl_output.type, WDL.Type.File):
                type_of = "File"
            # check for output with referee and referee expr of type WDL.Expr.Apply
            if not wdl_output.expr:
                raise ValueError("Missing expression")
            get_expr = cast(WDL.Expr.Get, wdl_output.expr)
            expr_ident = cast(WDL.Expr.Ident, get_expr.expr)
            expr_ident_name = expr_ident.name
            if expr_ident.referee and isinstance(
                expr_ident.referee.expr, WDL.Expr.Apply
            ):
                reference_expr = expr_ident.referee
                ref_expr_to_apply = reference_expr.expr
                ref_function = ref_expr_to_apply.function_name
                ref_arguments = ref_expr_to_apply.arguments
                if ref_function == "_add":
                    # return true value for a javascript tenary
                    first_arg, second_arg = ref_arguments
                    second_arg_value = self.get_wdl_literal(second_arg.literal)
                    first_arg_fun_name = first_arg.function_name
                    if first_arg_fun_name == "basename":
                        only_argument = first_arg.arguments[0]
                        only_argument_expr_name = only_argument.expr.name
                    true_tenary = f"inputs.{only_argument_expr_name}.{first_arg_fun_name} + '{second_arg_value}'"

                glob_expr = (
                    "$(inputs."
                    + reference_expr.name
                    + " === null ?"
                    + f" ({true_tenary}) : inputs.{reference_expr.name})"
                )
            else:
                glob_expr = f"$(inputs.{expr_ident_name})"

            outputs.append(
                cwl.CommandOutputParameter(
                    id=output_name,
                    type=type_of,
                    outputBinding=cwl.CommandOutputBinding(glob=glob_expr),
                )
            )
        return outputs


def main() -> None:
    """Entry point."""
    # Command-line parsing.
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow", help="Path to WDL workflow")
    parser.add_argument("-o", "--output", help="Name of resultant CWL file")
    args = parser.parse_args()

    # write to a file in oop_cwl_files
    if args.output:
        with open(args.output, "w") as result:
            result.write(str(Converter.load_wdl_tree(args.workflow)))

    # Converter.load_wdl_tree("wdl2cwl/tests/wdl_files/bowtie_1.wdl")
    # Converter.load_wdl_tree("wdl2cwl/tests/wdl_files/bcftools_stats.wdl")


if __name__ == "__main__":

    main()
