cwlVersion: v1.2
id: Sample
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

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
hints:
  - class: ResourceRequirement
    outdirMin: 1024
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
baseCommand:
  - bash
  - script.bash
outputs:
  - id: subsampledReads
    type: File
    outputBinding:
        glob: $(inputs.outFilePath)
