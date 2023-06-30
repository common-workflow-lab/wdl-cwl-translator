cwlVersion: v1.2
id: Echo
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            echo "Hello world"
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
  - class: ResourceRequirement
    coresMin: 8
    outdirMin: 1024
inputs: []
baseCommand:
  - bash
  - script.bash
outputs:
  - id: out
    type: stdout
