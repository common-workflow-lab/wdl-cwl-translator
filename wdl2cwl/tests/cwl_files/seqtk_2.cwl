class: CommandLineTool
id: Sample
inputs:
  - id: sequenceFile
    type: File
  - id: fractionOrNumber
    type: float
  - id: zip
    type:
      - boolean
      - 'null'
  - id: preCommand
    type:
      - string
      - 'null'
  - id: seed
    type:
      - int
      - 'null'
  - id: outFilePath
    default: subsampledReads.fq.gz
    type: string
  - id: twoPassMode
    default: false
    type: boolean
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
                    $(inputs["preCommand"] === null ? "" : inputs["preCommand"])
                    seqtk sample \
            	    $(inputs["seed"] === null ? "" : "-s " + inputs["seed"] ) \
                    $(inputs["twoPassMode"] ? "-2 " : "") \
                    $(inputs.sequenceFile.path) \
                    $(inputs.fractionOrNumber) \
                    $(inputs["zip"] === null ? "" : "| gzip") \
                    >  $(inputs.outFilePath)
  - class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
