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
  in:
  - id: input_file
  out:
  - id: output_file
  run: dockstore-tool-md5sum.cwl

