class: CommandLineTool
id: Echo
inputs: []
outputs:
  - id: out
    type: stdout
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

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
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
