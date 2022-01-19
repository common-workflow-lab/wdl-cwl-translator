class: CommandLineTool
id: Sample
inputs:
  - id: sequenceFile
    type: File
  - id: outFilePath
    default: subsampledReads.fq.gz
    type: string
  - id: twoPassMode
    default: false
    type: boolean
  - id: fractionOrNumber
    type: float
  - id: zip
    default: true
    type: boolean
  - id: preCommand
    type:
      - string
      - 'null'
  - id: seed
    type:
      - int
      - 'null'
outputs:
  - id: subsampledReads
    type: File
    outputBinding:
        glob: $(inputs.outFilePath)
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e -o pipefail
            mkdir -p "\$(dirname $(inputs.outFilePath))"
            $(inputs.preCommand)
            seqtk sample \
            $(inputs.seed === null ? "" : "-s " + inputs.seed) \
            $(inputs.twoPassMode ? "-2 " : "") \
            $(inputs.sequenceFile.path) \
            $(inputs.fractionOrNumber) \
            $(inputs.zip ? "| gzip" : "") \
            >  $(inputs.outFilePath)
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    outdirMin: 1024
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
