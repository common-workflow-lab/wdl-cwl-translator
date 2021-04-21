#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
#import cwl_utils.parser_v1_2 as cwl

from cwl_utils import parser_v1_2 as cwl
from ruamel import yaml

def main() -> None:
    """Generate a CWL object to match "cat-tool.cwl"."""
    inputs = [cwl.WorkflowInputParameter(id="input_file", type="File",)]
    outputs = [
        cwl.WorkflowOutputParameter(
            id="output_file",
            type="File",
            outputSource='md5sum/output_file',
            secondaryFiles='md5sum.cwl'
        )
    ]

    stepsIn = [
        cwl.WorkflowStepInput(
            id="input_file",
        )
    ]

    stepsOut = [
        cwl.WorkflowStepOutput(
            id="output_file",
        )
    ]
    steps = [
        cwl.WorkflowStep(
            id = "md5sum",
            run = "dockstore-tool-md5sum.cwl",
            in_ = stepsIn,
            out = stepsOut,
        )
    ]

    cat_tool = cwl.Workflow(
        inputs = inputs,
        outputs = outputs,
        steps = steps,
        cwlVersion= "v1.0",
    )

    print(yaml.main.round_trip_dump(cat_tool.save()))
 
if __name__ == "__main__":
    main()

'''python create_cwl_from_objects.py > result.cwl
cwltool --validate result.cwl'''
