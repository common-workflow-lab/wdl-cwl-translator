class: CommandLineTool
id: Format
inputs:
  - id: inputFiles
    type:
      - items: File
        type: array
  - id: format
    default: fasta
    type: string
  - id: outputPath
    default: seq_data.sdf
    type: string
  - id: rtgMem
    default: 8G
    type: string
  - id: memory
    default: 9G
    type: string
  - id: dockerImage
    default: quay.io/biocontainers/rtg-tools:3.10.1--0
    type: string
outputs:
  - id: sdf
    type: File
    outputBinding:
        glob: outputPath
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p \$(dirname $(inputs.outputPath))
            rtg RTG_MEM=$(inputs.rtgMem) format -f $(inputs.format) \
            -o $(inputs.outputPath) \
  - class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
