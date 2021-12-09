import os
from typing import List, Union, Optional, Callable
import WDL
import cwl_utils.parser.cwl_v1_2 as cwl

from io import StringIO
import textwrap
import re

from ruamel.yaml import scalarstring
from ruamel.yaml.main import YAML



class Converter:

    @staticmethod
    def load_wdl_tree(doc: str):
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
        if isinstance(obj, WDL.Task):
            return self.load_wdl_task(obj)
        elif isinstance(obj, WDL.Workflow):
            return self.load_wdl_workflow(obj)

    def load_wdl_workflow(self, obj: WDL.Workflow):
        print(f"Workflow {obj.name} loaded")
    
    def load_wdl_task(self, obj: WDL.Task):
        runtime = obj.runtime
        
        cwl_inputs = self.get_cwl_inputs(obj.inputs)
        base_command = ["bash", "example.sh"]

        cat_tool = cwl.CommandLineTool(
        id=obj.name,
        inputs=cwl_inputs,
        requirements=None,
        outputs=[],
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

    def get_cwl_inputs(self, wdl_inputs: List[str]):
        inputs = []

        for wdl_input in wdl_inputs:
            input_name = wdl_input.name
            input_value = None

            if isinstance(wdl_input, WDL.Type.Array):
                input_type = 'File'
                type_of = [cwl.CommandInputArraySchema(items=input_type, type="array")]
            elif isinstance(wdl_input, WDL.Type.String):
                type_of = "string"
            elif isinstance(wdl_input, WDL.Type.Boolean):
                type_of = "boolean"
            elif isinstance(wdl_input, WDL.Type.Int):
                type_of = "int"
            else:
                print(wdl_input.type)
                type_of = ""


            if wdl_input.type.optional: type_of = [type_of, "null"]

            if wdl_input.expr is not None: input_value = wdl_input.expr.literal.value

            inputs.append(cwl.CommandInputParameter(id=input_name, type=type_of, default=input_value))

        return inputs
        

    def translate_command(self, expr: WDL.Expr.Base, inputs: List[str]):

        if expr is None:
            return None
        
        if isinstance(expr, WDL.Expr.Array):
            return [self.translate_expr(e) for e in expr.items]

        if isinstance(expr, WDL.Expr.String):
            return self.translate_command_string(expr)
        elif isinstance(expr, (WDL.Expr.Int, WDL.Expr.Boolean, WDL.Expr.Float)):
            return self.literal.value
        if isinstance(expr, WDL.Expr.Placeholder):
            return self.translate_expr(expr.expr)



    def translate_command_string(self, string: WDL.Expr.String):
        # print("this is the literal", string.literal)
        # if string.literal is not None:
        #     return str(string.literal).lstrip('"').rstrip('"')

        # elements = {}
        # counter = 1
        # _format = str(string).lstrip('"').rstrip('"')
        # print("from _format", _format)

        # for placeholder in string.children:


        #     print(placeholder)
        pass




if __name__ == '__main__':
    import sys
    import argparse
    
    
    # Command-line parsing.
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow", help="Path to WDL workflow")

    # args = parser.parse_args()


    try:
        converted = Converter.load_wdl_tree("wdl2cwl/tests/wdl_files/bowtie_1.wdl")
        # converted = Converter.load_wdl_tree(args.workflow)
    except WDL.Error.SyntaxError as err:
        print(err)
    except WDL.Error.ValidationError as err:
        print(err)
    except WDL.Error.MultipleValidationErrors as err:
        for error in err.exceptions:
            print(error)
