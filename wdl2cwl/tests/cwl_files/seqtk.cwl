class: CommandLineTool
id: Sample
inputs:
  - id: sequenceFile
    type: File
  - id: fractionOrNumber
    type: float
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
  - id: zip
    default: true
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
                    mkdir -p "\$(dirname $(inputs["outFilePath"]===null?"":inputs["outFilePath"]["class"]==="File"? "$(inputs.outFilePath.path)": "$(inputs.outFilePath)"))"
                    $(inputs["preCommand"]===null?"":inputs["preCommand"]["class"]==="File"? "$(inputs.preCommand.path)": "$(inputs.preCommand)")
                    seqtk sample \
            	    -s $(inputs.seed) \
                    $(inputs["sequenceFile"]===null?"":inputs["sequenceFile"]["class"]==="File"? "$(inputs.sequenceFile.path)": "$(inputs.sequenceFile)") \
                    $(inputs["fractionOrNumber"]===null?"":inputs["fractionOrNumber"]["class"]==="File"? "$(inputs.fractionOrNumber.path)": "$(inputs.fractionOrNumber)") \
                    >  $(inputs["outFilePath"]===null?"":inputs["outFilePath"]["class"]==="File"? "$(inputs.outFilePath.path)": "$(inputs.outFilePath)")
  - class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
