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
  default: 'false'
  type: boolean
- id: zip
  default: 'true'
  type: boolean
outputs:
- id: subsampledReads
  type: File
  outputBinding:
    glob: outFilePath
requirements:
- class: InitialWorkDirRequirement
  listing:
  - entryname: example.sh
    entry: "set -e -o pipefail mkdir -p \"$(dirname $(inputs.outFilePath))\" $(inputs.preCommand)\
      \ seqtk sample \\-s $(inputs.seed) \\$(if ['${inputs.twoPassMode}' eq 'true'\
      \ ] then echo-2  else echo fi ) \\$(inputs.sequenceFile.path) \\$(inputs.fractionOrNumber)\
      \ \\$(if ['${inputs.zip}' eq 'true' ] then echo| gzip else echo fi ) \\>  $(inputs.outFilePath) "
cwlVersion: v1.2
baseCommand:
- sh
- example.sh
