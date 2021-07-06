class: CommandLineTool
id: Sample
inputs:
  - id: sequenceFile
    type: File
  - id: fractionOrNumber
    type: float
  - id: outFilePath
    default: subsampledReads.fq.gz
    type: string
  - id: twoPassMode
    default: false
    type: boolean
  - id: zip
    default: true
    type: boolean
outputs:
  - id: quality_yield_metrics
    type: File
    outputBinding:
        glob: '"abc"'
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e -o pipefail
            mkdir -p "$(dirname $(inputs.outFilePath))"
            $(inputs.preCommand)
            seqtk sample \
            $(inputs."-s "+seed) \
            $(inputs.true="-2 "false=""twoPassMode) \
            $(inputs.sequenceFile.path) \
            $(inputs.fractionOrNumber) \
            $(inputs.true="| gzip"false=""zip) \
            >  $(inputs.outFilePath)
  - class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
