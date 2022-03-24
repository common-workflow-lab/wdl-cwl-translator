"""Main entrypoint for WDL2CWL."""
import argparse
import re
import textwrap
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)

import cwl_utils.parser.cwl_v1_2 as cwl
import regex  # type: ignore
import sys
from ruamel.yaml import scalarstring
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml.main import YAML

import WDL
import WDL.CLI
import WDL._parser  # delete when reloading bug is fixed upstream
from wdl2cwl import _logger
from wdl2cwl.errors import WDLSourceLine, ConversionException
from wdl2cwl.expr import (
    get_expr,
    get_literal_value,
    get_step_input_expr,
    get_expr_name,
    translate_wdl_placeholder,
)
from wdl2cwl.util import get_mem_in_bytes, get_input, ConversionContext


def convert(doc: str) -> Dict[str, Any]:
    """Convert a WDL workflow, reading the file, into a CWL workflow Python object."""
    WDL._parser._lark_comments_buffer.clear()
    root_folder = str(Path(doc).parent)
    try:
        doc_tree = WDL.load(
            doc, [root_folder], read_source=WDL.CLI.make_read_source(False), check_quant=True  # type: ignore[no-untyped-call]
        )
    except (
        WDL.Error.SyntaxError,
        WDL.Error.ImportError,
        WDL.Error.ValidationError,
        WDL.Error.MultipleValidationErrors,
    ) as exn:
        WDL.CLI.print_error(exn)
        raise exn

    parser = Converter()
    if doc_tree.workflow:
        return sort_cwl(parser.load_wdl_workflow(doc_tree.workflow).save())
    if len(doc_tree.tasks) == 1:
        return sort_cwl(parser.load_wdl_objects(doc_tree.tasks[0]).save())
    else:
        return {
            "cwlVersion": "v1.2",
            "$graph": [
                sort_cwl(parser.load_wdl_objects(task).save())
                for task in doc_tree.tasks
            ],
        }


def sort_cwl(document: Dict[str, Any]) -> CommentedMap:
    """Sort the sections of the CWL document in a more meaningful order."""
    keyorder = [
        "cwlVersion",
        "id",
        "class",
        "label",
        "doc",
        "requirements",
        "hints",
        "inputs",
        "stdin",
        "baseCommand",
        "steps",
        "expression",
        "arguments",
        "stderr",
        "stdout",
        "outputs",
        "successCodes",
        "temporaryFailCodes",
        "permanentFailCodes",
    ]
    return CommentedMap(
        sorted(
            document.items(),
            key=lambda i: keyorder.index(i[0]) if i[0] in keyorder else 100,
        )
    )


CWLArrayTypes = TypeVar(
    "CWLArrayTypes",
    cwl.InputArraySchema,
    cwl.OutputArraySchema,
    cwl.CommandInputArraySchema,
    cwl.CommandOutputArraySchema,
)
CWLRecordTypes = TypeVar(
    "CWLRecordTypes",
    cwl.InputRecordSchema,
    cwl.OutputRecordSchema,
    cwl.CommandInputRecordSchema,
    cwl.CommandOutputRecordSchema,
)


def get_cwl_docker_requirements(
    wdl_docker: Union[WDL.Expr.Get, WDL.Expr.String]
) -> Optional[cwl.ProcessRequirement]:
    """Translate WDL Runtime Docker requirements to CWL Docker Requirement."""
    if isinstance(wdl_docker, WDL.Expr.String) and wdl_docker.literal:
        dockerpull = get_literal_value(wdl_docker)
    else:
        if isinstance(wdl_docker, WDL.Expr.String):
            parts = wdl_docker.parts
            docker_placeholder = [
                pl_holder
                for pl_holder in parts
                if isinstance(pl_holder, WDL.Expr.Placeholder)
            ]
            dockerpull_expr = docker_placeholder[0].expr.expr  # type: ignore[attr-defined]
        else:
            dockerpull_expr = wdl_docker.expr
        if dockerpull_expr is None or not isinstance(dockerpull_expr, WDL.Expr.Ident):
            raise WDLSourceLine(wdl_docker, ConversionException).makeError(
                f"Unsupported type: {type(dockerpull_expr)}: {dockerpull_expr}"
            )
        dockerpull_referee = dockerpull_expr.referee
        if dockerpull_referee is None:
            raise WDLSourceLine(wdl_docker, ConversionException).makeError(
                f"Unsupported type: {type(dockerpull_referee)}"
            )
        if dockerpull_referee.expr:
            dockerpull = get_literal_value(dockerpull_referee.expr)
        else:
            _logger.warning(f"Unable to extract docker reference from {wdl_docker}")
            return None
    return cwl.DockerRequirement(dockerPull=dockerpull)


read_funcs = {
    "read_string": r"$(self[0].contents.replace(/[\r\n]+$/, ''))",
    "read_int": "$(parseInt(self[0].contents))",
    "read_float": "$(parseFloat(self[0].contents))",
    "read_boolean": """${
  var contents = self[0].contents.trim().toLowerCase()
  if (contents == 'true') { return true;}
  if (contents == 'false') { return false;}
  throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
}""",
    "read_lines": r"""${
  var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
  // ^ remove any trailing newline to prevent a null being returned
  return contents.split(/\r\n|\r|\n/);
}""",
    "read_tsv": r"""${
  var result = Array();
  var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
  // ^ remove any trailing newline to prevent a null being returned
  contents.split(/\r\n|\r|\n/).forEach(function(line) {
    result.push(line.split('\t'));
  });
  return result;
}""",
}


class Converter:
    """Object that handles WDL Workflows and task conversion to CWL."""

    def __init__(self) -> None:
        """Initialize the context used by the object to prevent inconsistent behaviours."""
        self.ctx = ConversionContext()

    def load_wdl_objects(
        self, obj: Union[WDL.Tree.Task, WDL.Tree.Workflow]
    ) -> Union[cwl.CommandLineTool, cwl.Workflow]:
        """Load a WDL SourceNode obj and returns either a Task or a Workflow."""
        if isinstance(obj, WDL.Tree.Task):
            return self.load_wdl_task(obj)
        elif isinstance(obj, WDL.Tree.Workflow):
            return self.load_wdl_workflow(obj)

    def load_wdl_workflow(self, obj: WDL.Tree.Workflow) -> cwl.Workflow:
        """Load WDL workflow and convert to CWL."""
        wf_name = obj.name
        requirements: List[cwl.ProcessRequirement] = [cwl.InlineJavascriptRequirement()]
        inputs = self.get_cwl_workflow_inputs(obj.available_inputs, obj.parameter_meta)
        outputs = [
            cwl.WorkflowOutputParameter(
                id=f"{wf_name}.{output_id}",
                type=output_type,
                outputSource=output_source,
                doc=obj.meta.get(output_id, None),
            )
            for output_id, output_type, output_source in self.get_workflow_outputs(
                obj.effective_outputs
            )
        ]
        wf_steps: List[cwl.WorkflowStep] = []
        wf_description = (
            obj.meta.pop("description") if "description" in obj.meta else None
        )
        scatter_present = False
        step_valuefrom = False
        for body_part in obj.body:
            if not isinstance(body_part, (WDL.Tree.Call, WDL.Tree.Scatter)):
                _logger.warning(
                    WDLSourceLine(body_part).makeError(
                        "Warning: unhandled Workflow node type:"
                    )
                    + " %s",
                    type(body_part),
                )
                continue
            with WDLSourceLine(body_part, ConversionException):
                if isinstance(body_part, WDL.Tree.Call):
                    step = self.get_workflow_call(body_part)
                    wf_steps.append(step)
                    for inp in step.in_:
                        if inp.valueFrom is not None:
                            step_valuefrom = True
                elif isinstance(body_part, WDL.Tree.Scatter):
                    scatter_present = True
                    wf_steps.extend(self.get_workflow_scatter(body_part))
        if scatter_present:
            requirements.append(cwl.ScatterFeatureRequirement())
        if step_valuefrom:
            requirements.append(cwl.StepInputExpressionRequirement())

        return cwl.Workflow(
            id=wf_name,
            cwlVersion="v1.2",
            doc=wf_description,
            requirements=requirements if requirements else None,
            inputs=inputs,
            steps=wf_steps,
            outputs=outputs,
        )

    def get_workflow_scatter(self, scatter: WDL.Tree.Scatter) -> List[cwl.WorkflowStep]:
        """Get the CWL Workflow Step equivalent of a list of WDL Scatter Object."""
        scatter_steps: List[cwl.WorkflowStep] = []
        for body_obj in scatter.body:
            with WDLSourceLine(body_obj, ConversionException):
                if isinstance(body_obj, WDL.Tree.Call):
                    wf_step = self.get_workflow_call(body_obj)
                wf_step.scatter = self.ctx.scatter_names.pop()
                scatter_steps.append(wf_step)
        return scatter_steps

    def get_workflow_call(self, call: WDL.Tree.Call) -> cwl.WorkflowStep:
        """Get the CWL Workflow Step equivalent of a WDL Call Object."""
        callee = call.callee
        if callee is None:
            raise ConversionException("callee is None. This should not be possible.")
        cwl_callee_inputs = self.get_cwl_task_inputs(callee.inputs)[0]
        call_inputs = call.inputs
        inputs_from_call: Dict[str, Tuple[str, Dict[str, Any]]] = {}
        input_defaults = set()
        if call_inputs:
            scatter_handled = False
            for key, value in call_inputs.items():
                with WDLSourceLine(value, ConversionException):
                    if not isinstance(value, (WDL.Expr.Get, WDL.Expr.Apply)):
                        input_defaults.add(key)
                    if (
                        isinstance(value, WDL.Expr.Apply)
                        and value.function_name == "select_all"
                        and isinstance(value.arguments[0], WDL.Expr.Array)
                        and isinstance(
                            value.arguments[0].items[0],
                            (WDL.Expr.Get, WDL.Expr.String),
                        )
                    ):
                        input_expr, value_from = get_step_input_expr(
                            value.arguments[0].items[0], self.ctx
                        )
                        value_str = f"$([{value_from if value_from else 'self'}])"
                        inputs_from_call[key] = (
                            input_expr.replace(".", "/"),
                            {
                                "pickValue": "all_non_null",
                                "valueFrom": value_str,
                            },
                        )
                    else:
                        input_expr, value_from = get_step_input_expr(value, self.ctx)  # type: ignore[arg-type]
                        inputs_from_call[key] = (
                            input_expr.replace(".", "/")
                            if isinstance(input_expr, str)
                            else input_expr,
                            {"valueFrom": f"$({value_from})"} if value_from else {},
                        )
                    if self.ctx.scatter_names and not scatter_handled:
                        self.ctx.scatter_names[-1] = key
                        scatter_handled = True
        wf_step_inputs: List[cwl.WorkflowStepInput] = []
        for inp in cwl_callee_inputs:
            call_inp_id = f"{call.name}.{inp.id}"
            source_str, extras = inputs_from_call.get(
                cast(str, inp.id), (call_inp_id, {})
            )

            if inp.id not in input_defaults:
                wf_step_inputs.append(
                    cwl.WorkflowStepInput(id=inp.id, source=source_str, **extras)
                )
            else:
                wf_step_inputs.append(
                    cwl.WorkflowStepInput(id=inp.id, default=source_str, **extras)
                )
        wf_step_outputs = (
            [cwl.WorkflowStepOutput(id=output.name) for output in callee.outputs]
            if callee.outputs
            else []
        )
        wf_step_run = self.load_wdl_objects(callee)
        wf_step = cwl.WorkflowStep(
            wf_step_inputs,
            id=call.name,
            run=wf_step_run,
            out=wf_step_outputs,
        )
        return wf_step

    def load_wdl_task(self, obj: WDL.Tree.Task) -> cwl.CommandLineTool:
        """Load task and convert to CWL."""
        description = obj.meta.pop("description", None)
        cwl_inputs, restrictions = self.get_cwl_task_inputs(
            obj.inputs, obj.parameter_meta
        )
        run_script = False
        hints, requirements = self.get_cwl_hints_and_requirements(obj)
        for req in requirements:
            if isinstance(req, cwl.InitialWorkDirRequirement):
                run_script = True
        arguments = [{"valueFrom": "${" + item + "}"} for item in restrictions]
        tool = cwl.CommandLineTool(
            id=obj.name,
            doc=description,
            inputs=cwl_inputs,
            outputs=None,
            hints=hints,
            requirements=requirements,
            cwlVersion="v1.2",
            baseCommand=["bash", "script.bash"] if run_script else ["true"],
            arguments=arguments if arguments else None,
        )
        tool.outputs = self.set_cwl_task_outputs(obj.outputs, obj.parameter_meta, tool)
        if obj.meta:
            _logger.warning("Skipping meta: %s", obj.meta)
        return tool

    def get_cwl_hints_and_requirements(
        self, obj: WDL.Tree.Task
    ) -> Tuple[List[cwl.ProcessRequirement], List[cwl.ProcessRequirement]]:
        """Produce the CWL Requirements list."""
        command = obj.command
        runtime = obj.runtime
        requirements: List[cwl.ProcessRequirement] = []
        hints: List[cwl.ProcessRequirement] = []
        if "docker" in runtime:
            with WDLSourceLine(runtime["docker"], ConversionException):
                docker = get_cwl_docker_requirements(
                    runtime["docker"]  # type: ignore[arg-type]
                )
                if docker:
                    hints.append(docker)
        command_req = self.get_cwl_command_requirements(command.parts)
        if command_req:
            requirements.append(command_req)
        requirements.append(cwl.InlineJavascriptRequirement())
        requirements.append(cwl.NetworkAccess(networkAccess=True))
        cpu_requirement = (
            self.get_cpu_requirement(runtime["cpu"]) if "cpu" in runtime else None
        )
        if "memory" in runtime:
            with WDLSourceLine(runtime["memory"], ConversionException):
                memory_requirement = self.get_memory_requirement(
                    runtime["memory"]  # type: ignore[arg-type]
                )
        else:
            memory_requirement = None
        if "disks" in runtime:
            with WDLSourceLine(runtime["memory"], ConversionException):
                outdir_requirement = self.get_outdir_requirement(
                    runtime["disks"]  # type: ignore[arg-type]
                )
        else:
            outdir_requirement = 1024
        hints.append(
            cwl.ResourceRequirement(
                coresMin=cpu_requirement,
                ramMin=memory_requirement,
                outdirMin=outdir_requirement,
            )
        )
        if "time_minutes" in runtime:
            with WDLSourceLine(runtime["time_minutes"], ConversionException):
                time_minutes = self.get_time_minutes_requirement(
                    runtime["time_minutes"]  # type: ignore[arg-type]
                )
            hints.append(
                cwl.ToolTimeLimit(
                    timelimit=time_minutes,
                )
            )
        return hints, requirements

    def get_time_minutes_requirement(
        self, time_minutes: WDL.Expr.Get
    ) -> Union[str, int]:
        """Produce the time limit expression from WDL runtime time minutes."""
        with WDLSourceLine(time_minutes, ConversionException):
            if isinstance(time_minutes, (WDL.Expr.Int, WDL.Expr.Float)):
                return cast(int, get_literal_value(time_minutes)) * 60
            time_minutes_str = get_expr(time_minutes, self.ctx)
        return f"$({time_minutes_str} * 60)"

    def get_memory_requirement(
        self, memory_runtime: Union[WDL.Expr.Get, WDL.Expr.String]
    ) -> Union[str, float]:
        """Translate WDL Runtime Memory requirement to CWL Resource Requirement."""
        with WDLSourceLine(memory_runtime, ConversionException):
            if isinstance(memory_runtime, WDL.Expr.String):
                ram_min_literal = self.get_memory_literal(memory_runtime)
                return ram_min_literal
            elif isinstance(memory_runtime, WDL.Expr.Apply):
                expr, unit = memory_runtime.arguments
                ram_min = get_expr(expr, self.ctx)
                return self.get_ram_min_js(
                    ram_min, str(get_literal_value(unit)).strip()
                )
            ram_min = get_expr_name(memory_runtime.expr)  # type: ignore[arg-type]
            return self.get_ram_min_js(ram_min, "")

    def get_memory_literal(self, memory_runtime: WDL.Expr.String) -> Union[float, str]:
        """Get the literal value for memory requirement with type WDL.Expr.String."""
        if memory_runtime.literal is None:
            if len(memory_runtime.parts) == 4:
                _, amount, unit, _ = memory_runtime.parts
            else:
                _, amount, _, unit, _ = memory_runtime.parts
            if isinstance(amount, WDL.Expr.Placeholder):
                amount_str = get_expr(amount, self.ctx)
            else:
                amount_str = amount
            if isinstance(unit, WDL.Expr.Placeholder):
                unit_str = get_expr(unit, self.ctx)
            else:
                unit_str = unit.strip()
            return self.get_ram_min_js(amount_str, unit_str)

        ram_min = str(get_literal_value(memory_runtime))
        unit_result = re.search(r"[a-zA-Z]+", ram_min)
        if not unit_result:
            raise ConversionException("Missing Memory units, yet still a string?")
        unit = unit_result.group()
        value = float(ram_min.split(unit)[0])
        unit_base, unit_exponent = get_mem_in_bytes(unit)
        return float((value * (unit_base**unit_exponent)) / (2**20))

    def get_outdir_requirement(
        self, outdir: Union[WDL.Expr.Get, WDL.Expr.Apply, WDL.Expr.Placeholder]
    ) -> Union[int, str]:
        """Produce the memory requirement for the output directory from WDL runtime disks."""
        with WDLSourceLine(outdir, ConversionException):
            if (
                isinstance(outdir, (WDL.Expr.Apply, WDL.Expr.String))
                and outdir.literal is None
            ):
                if (
                    isinstance(outdir, WDL.Expr.Apply)
                    and not outdir.function_name == "_add"
                ):
                    # If it contains an apply expr we don't want to process the _add function
                    # that concatenates it to the chars in the string
                    expr_str = get_expr(outdir, self.ctx)
                    if isinstance(outdir.type, WDL.Type.String):
                        return f"$(parseFloat({expr_str}) * 1024)"
                    else:
                        return f"$(({expr_str}) * 1024)"
                # apply exprs contain arguments and strings contain parts both are lists
                # they could contain expressions that represent the runtime disk
                list_object = (
                    outdir.arguments
                    if hasattr(outdir, "arguments")
                    else outdir.parts  # type: ignore[attr-defined]
                )
                for obj in list_object:
                    if isinstance(
                        obj, (WDL.Expr.Get, WDL.Expr.Apply, WDL.Expr.Placeholder)
                    ):
                        # avoid python strings only WDL expressions are handled.
                        return self.get_outdir_requirement(obj)
            elif isinstance(outdir, (WDL.Expr.Get, WDL.Expr.Placeholder)):
                expr_str = get_expr(outdir, self.ctx)
                return (
                    f"$(({expr_str}) * 1024)"
                    if not expr_str.isdigit()
                    else int(expr_str) * 1024
                )
            value = re.search(r"[0-9]+", str(get_literal_value(outdir))).group()  # type: ignore[union-attr]

            return int(value) * 1024

    def get_ram_min_js(self, ram_min_ref_name: str, unit: str) -> str:
        """Get memory requirement for user input."""
        append_str: str = ""
        if unit:
            if "inputs." in unit:
                append_str = "${\nvar unit = " + unit + ";"
            else:
                append_str = '${\nvar unit = "' + unit + '";'
        else:
            append_str = (
                "${\nvar unit = " + ram_min_ref_name + '.match(/[a-zA-Z]+/g).join("");'
            )
        js_str = (
            append_str
            + "\nvar value = parseInt(`${"
            + ram_min_ref_name
            + "}`.match(/[0-9]+/g));\n"
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
            + 'else throw "Unknown units: " + unit;\n'
            + "return parseInt(memory);\n}"
        )

        return js_str

    def get_cwl_type(
        self,
        input_type: WDL.Type.Base,
        array_type: Type[CWLArrayTypes],
        record_type: Type[CWLRecordTypes],
    ) -> Union[str, CWLArrayTypes, CWLRecordTypes]:
        """Determine the CWL type for a WDL input declaration."""
        if isinstance(input_type, WDL.Type.File):
            return "File"
        elif isinstance(input_type, WDL.Type.String):
            return "string"
        elif isinstance(input_type, WDL.Type.Boolean):
            return "boolean"
        elif isinstance(input_type, WDL.Type.Int):
            return "int"
        elif isinstance(input_type, WDL.Type.Float):
            return "float"
        elif isinstance(input_type, WDL.Type.Array):
            sub_type = self.get_cwl_type(input_type.item_type, array_type, record_type)
            return array_type(type="array", items=sub_type)
        elif isinstance(input_type, WDL.Type.StructInstance):
            return record_type(
                type="record",
                name=input_type.type_name,
                fields=self.get_struct_inputs(input_type.members),
            )
        else:
            raise WDLSourceLine(input_type, ConversionException).makeError(
                f"Input of type {input_type} is not yet handled."
            )

    def get_workflow_outputs(
        self,
        outputs: WDL.Env.Bindings[WDL.Type.Base],
    ) -> Iterator[
        Tuple[str, Union[cwl.OutputArraySchema, cwl.OutputRecordSchema, str], str]
    ]:
        """Return the name, CWL type, and source for a workflow's effective_outputs()."""
        for item in outputs:
            with WDLSourceLine(item.info, ConversionException):
                output_name = item.name
                item_expr = item.info.expr
                output_source = item_expr.expr.name[::-1].replace(".", "/", 1)[::-1]
                # replace just the last occurrence of a period with a slash
                # by first reversing the string and the replace the first occurrence
                # then reversing the result
                wdl_output = item.info
                type_of = self.get_cwl_type(
                    wdl_output.type, cwl.OutputArraySchema, cwl.OutputRecordSchema
                )
                yield (output_name, type_of, output_source)

    def get_cpu_requirement(self, cpu_runtime: WDL.Expr.Base) -> Union[int, float, str]:
        """Translate WDL Runtime CPU requirement to CWL Resource Requirement."""
        if isinstance(cpu_runtime, (WDL.Expr.Int, WDL.Expr.Float)):
            return cast(Union[int, float], get_literal_value(cpu_runtime))
        elif isinstance(cpu_runtime, WDL.Expr.String):
            literal_str = cast(Any, get_literal_value(cpu_runtime))
            return int(literal_str) if "." not in literal_str else float(literal_str)
        return f"$({get_expr(cpu_runtime, self.ctx)})"

    def get_cwl_command_requirements(
        self, wdl_commands: List[Union[str, WDL.Expr.Placeholder]]
    ) -> Optional[cwl.InitialWorkDirRequirement]:
        """Translate WDL commands into CWL Initial WorkDir Requirement."""
        command_str: str = ""
        for wdl_command in wdl_commands:
            if isinstance(wdl_command, str):
                command_str += wdl_command.replace("$(", "\\$(")
            else:
                command_str += f"$({translate_wdl_placeholder(wdl_command, self.ctx)})"
        command_str = textwrap.dedent(command_str)
        return (
            cwl.InitialWorkDirRequirement(
                listing=[cwl.Dirent(entry=command_str, entryname="script.bash")]
            )
            if len(command_str.strip()) > 0
            else None
        )

    def get_cwl_workflow_inputs(
        self,
        wdl_inputs: WDL.Env.Bindings[WDL.Tree.Decl],
        meta: Optional[Dict[str, Any]] = None,
    ) -> List[cwl.WorkflowInputParameter]:
        """Convert WDL inputs into CWL inputs and return a list of CWL Workflow Input Parameters."""
        inputs: List[cwl.WorkflowInputParameter] = []

        for input_decl in wdl_inputs:
            input_name = input_decl.name
            self.ctx.non_static_values.add(input_name)
            input_value = None
            wdl_input = input_decl.value
            type_of = self.get_cwl_type(
                wdl_input.type,
                cwl.InputArraySchema,
                cwl.InputRecordSchema,
            )

            if wdl_input.type.optional or isinstance(wdl_input.expr, WDL.Expr.Apply):
                final_type_of: Union[
                    List[
                        Union[
                            str,
                            cwl.InputArraySchema,
                            cwl.InputRecordSchema,
                        ]
                    ],
                    str,
                    cwl.InputArraySchema,
                    cwl.InputRecordSchema,
                ] = [type_of, "null"]
            else:
                final_type_of = type_of

            if wdl_input.expr is not None:
                with WDLSourceLine(wdl_input.expr, ConversionException):
                    input_value = get_literal_value(wdl_input.expr)
                    if input_value and final_type_of == "float":
                        input_value = float(input_value)

            doc: Optional[str] = None
            if meta and input_name in meta:
                if isinstance(meta[input_name], str):
                    doc = meta[input_name]
                elif "description" in meta[input_name]:
                    doc = meta[input_name]["description"]
            inputs.append(
                cwl.WorkflowInputParameter(
                    id=input_name,
                    type=final_type_of,
                    default=input_value,
                    doc=doc,
                )
            )

        return inputs

    def get_cwl_task_inputs(
        self,
        wdl_inputs: Optional[List[WDL.Tree.Decl]],
        meta: Optional[Dict[str, Any]] = None,
    ) -> Tuple[List[cwl.CommandInputParameter], List[str]]:
        """
        Convert WDL inputs into CWL inputs.

        Return a tuple: list of CWL Command Input Parameters and list of restriction checks.
        """
        inputs: List[cwl.CommandInputParameter] = []
        restriction_checks: List[str] = []

        if not wdl_inputs:
            return inputs, restriction_checks

        for wdl_input in wdl_inputs:
            input_name = wdl_input.name
            self.ctx.non_static_values.add(input_name)
            input_value = None

            type_of = self.get_cwl_type(
                wdl_input.type,
                cwl.CommandInputArraySchema,
                cwl.CommandInputRecordSchema,
            )

            if isinstance(wdl_input.type, WDL.Type.Array) and wdl_input.type.nonempty:
                restriction_checks.append(
                    f"if ({get_input(input_name)}.length == 0) "
                    "{throw "
                    f'"{input_name} must contain at least one item.";'
                    '} else { return "";}'
                )

            if wdl_input.type.optional or isinstance(wdl_input.expr, WDL.Expr.Apply):
                final_type_of: Union[
                    List[
                        Union[
                            str,
                            cwl.CommandInputArraySchema,
                            cwl.CommandInputRecordSchema,
                        ]
                    ],
                    str,
                    cwl.CommandInputArraySchema,
                    cwl.CommandInputRecordSchema,
                ] = [type_of, "null"]
                if isinstance(wdl_input.expr, WDL.Expr.Apply):
                    self.ctx.optional_cwl_null.add(input_name)
            else:
                final_type_of = type_of

            if wdl_input.expr is not None:
                if not isinstance(wdl_input.expr, WDL.Expr.Apply):
                    with WDLSourceLine(wdl_input.expr, ConversionException):
                        input_value = get_literal_value(wdl_input.expr)
                        if input_value and final_type_of == "float":
                            input_value = float(input_value)

            doc: Optional[str] = None
            if meta and input_name in meta:
                if isinstance(meta[input_name], str):
                    doc = meta[input_name]
                elif "description" in meta[input_name]:
                    doc = meta[input_name]["description"]
            inputs.append(
                cwl.CommandInputParameter(
                    id=input_name, type=final_type_of, default=input_value, doc=doc
                )
            )

        return inputs, restriction_checks

    def get_struct_inputs(
        self, members: Optional[Dict[str, WDL.Type.Base]]
    ) -> List[cwl.CommandInputRecordField]:
        """Get member items of a WDL struct and return a list of cwl.CommandInputRecordField."""
        inputs: List[cwl.CommandInputRecordField] = []
        if not members:
            return inputs
        for member, value in members.items():
            input_name = member
            type_of = self.get_cwl_type(
                value, cwl.CommandInputArraySchema, cwl.CommandInputRecordSchema
            )
            inputs.append(cwl.CommandInputRecordField(name=input_name, type=type_of))
        return inputs

    def set_cwl_task_outputs(
        self,
        wdl_outputs: List[WDL.Tree.Decl],
        meta: Optional[Dict[str, Any]],
        tool: cwl.CommandLineTool,
    ) -> List[cwl.CommandOutputParameter]:
        """Convert WDL outputs into CWL outputs and return a list of CWL Command Output Parameters."""
        outputs: List[cwl.CommandOutputParameter] = []

        if not wdl_outputs:
            return outputs

        for wdl_output in wdl_outputs:
            output_name = wdl_output.name
            doc: Optional[str] = None
            if meta and output_name in meta:
                if isinstance(meta[output_name], str):
                    doc = meta[output_name]
                else:
                    if "description" in meta[output_name]:
                        doc = meta[output_name]["description"]
            type_of = self.get_cwl_type(
                wdl_output.type,
                cwl.CommandOutputArraySchema,
                cwl.CommandOutputRecordSchema,
            )
            glob = type_of == "File" or (
                isinstance(type_of, cwl.CommandOutputArraySchema)
                and type_of.items == "File"
            )

            if not wdl_output.expr:
                raise WDLSourceLine(wdl_output, ConversionException).makeError(
                    "Missing expression"
                )

            if (
                isinstance(wdl_output.expr, WDL.Expr.Apply)
                and wdl_output.expr.function_name in read_funcs
            ):
                if (
                    isinstance(wdl_output.expr.arguments[0], WDL.Expr.Apply)
                    and wdl_output.expr.arguments[0].function_name == "stdout"
                ):
                    glob_str = tool.stdout = "_stdout"
                else:
                    glob_expr = get_expr(wdl_output.expr, self.ctx)
                    is_literal = wdl_output.expr.arguments[0].literal
                    if is_literal:
                        glob_str = str(get_literal_value(wdl_output.expr.arguments[0]))
                    else:
                        glob_str = f"$({glob_expr})"

                outputs.append(
                    cwl.CommandOutputParameter(
                        id=output_name,
                        doc=doc,
                        type=type_of,
                        outputBinding=cwl.CommandOutputBinding(
                            glob=glob_str,
                            loadContents=True,
                            outputEval=read_funcs[wdl_output.expr.function_name],
                        ),
                    )
                )
            elif (
                isinstance(wdl_output.expr, WDL.Expr.Apply)
                and wdl_output.expr.function_name == "stdout"
            ):
                outputs.append(
                    cwl.CommandOutputParameter(
                        id=output_name,
                        doc=doc,
                        type="stdout",
                    ),
                )
            else:
                if wdl_output.type.optional:
                    final_type_of: Union[
                        List[
                            Union[
                                str,
                                cwl.CommandOutputArraySchema,
                                cwl.CommandOutputRecordSchema,
                            ]
                        ],
                        str,
                        cwl.CommandInputArraySchema,
                        cwl.CommandOutputArraySchema,
                        cwl.CommandOutputRecordSchema,
                    ] = [type_of, "null"]
                else:
                    final_type_of = type_of
                if glob:
                    globs: List[str] = []
                    if isinstance(wdl_output.expr, WDL.Expr.Array):
                        targets: Sequence[WDL.Expr.Base] = wdl_output.expr.items
                    else:
                        targets = [wdl_output.expr]
                    for entry in targets:
                        glob_str = f"$({get_expr(entry, self.ctx)})"
                        if (
                            isinstance(wdl_output.expr, WDL.Expr.String)
                            and wdl_output.expr.literal is not None
                        ):
                            glob_str = glob_str[3:-2]
                        globs.append(glob_str)
                    final_glob = globs if len(globs) > 1 else globs[0]
                    outputs.append(
                        cwl.CommandOutputParameter(
                            id=output_name,
                            doc=doc,
                            type=final_type_of,
                            outputBinding=cwl.CommandOutputBinding(glob=final_glob),
                        )
                    )
                else:
                    outputEval = f"$({get_expr(wdl_output.expr, self.ctx)})"
                    outputs.append(
                        cwl.CommandOutputParameter(
                            id=output_name,
                            doc=doc,
                            type=final_type_of,
                            outputBinding=cwl.CommandOutputBinding(
                                outputEval=outputEval
                            ),
                        )
                    )
        return outputs


def main(args: Union[List[str], None] = None) -> None:
    """Entry point."""
    # Command-line parsing.
    parser = argparse.ArgumentParser(
        description="Converts WDL workflows into CWL workflows. Outputs "
        "to <stdout> by default."
    )
    parser.add_argument("workflow", help="Path to WDL workflow")
    parser.add_argument("-o", "--output", help="Name of output CWL file")
    parsed_args = parser.parse_args(args)

    cwl_result = convert(parsed_args.workflow)

    # Serialize result in YAML to either <stdout> or specified output file.
    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent = 4
    yaml.block_seq_indent = 2
    scalarstring.walk_tree(cwl_result)

    if parsed_args.output is None:
        yaml.dump(cwl_result, sys.stdout)
    else:
        with open(parsed_args.output, "w") as f:
            yaml.dump(cwl_result, f)


if __name__ == "__main__":

    main(sys.argv[1:])


__all__ = ["main", "Converter", "ConversionException"]
