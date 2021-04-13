#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import cwl_utils.parser_v1_2 as cwl

from ruamel import yaml

def main() -> None:
    """Generate a CWL object to match "cat-tool.cwl"."""
    inputs = [cwl.WorkflowInputParameter(id="input_file", type="File",)]
    outputs = [
        cwl.WorkflowOutputParameter(
            id="output_file",
            type="File",
            outputSource="md5sum/output_file",
            secondaryFiles='md5sum.cwl'
        )
    ]
    steps = [
        cwl.WorkflowStep(
            id = "md5sum",
            run = "dockstore-tool-md5sum.cwl",
            in_ = "input_file: input_file",
            out = "[output_file]"
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

'''
Expected Output:

cwlVersion: v1.0
class: Workflow

inputs:
  input_file: File

outputs:
  output_file:
    type: File
    outputSource: md5sum/output_file

steps:
  md5sum:
    run: dockstore-tool-md5sum.cwl
    in:
      input_file: input_file
    out: [output_file]

Actual Output:

class: Workflow
inputs:
- id: input_file
  type: File
outputs:
- id: output_file
  secondaryFiles: md5sum.cwl
  outputSource: md5sum/output_file
  type: File
cwlVersion: v1.0
steps:
- id: md5sum
  in: 'input_file: input_file'
  out: '[output_file]'
  run: dockstore-tool-md5sum.cwl

#
'''   
