class: CommandLineTool
id: Sample
inputs:
  - id: sequenceFile
    type: File
  - id: fractionOrNumber
    type: float
  - id: preCommand
    type: string
  - id: seed
    type: int
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
                    mkdir -p "\$(dirname $(inputs.outFilePath))"
                    $(inputs.preCommand)
                    seqtk sample \
            	    -s $(inputs.seed) \
                    ${var value="";
            if(inputs["twoPassMode"]){
            value="-2 ";
            }
            else{
            value="";
            }
            return value;} \
                    $(inputs.sequenceFile.path) \
                    $(inputs.fractionOrNumber) \
                    ${var value="";
            if(inputs["zip"]){
            value="| gzip";
            }
            else{
            value="";
            }
            return value;} \
                    >  $(inputs.outFilePath)
  - class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
