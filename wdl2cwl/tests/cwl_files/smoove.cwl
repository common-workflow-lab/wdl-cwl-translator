class: CommandLineTool
id: Call
inputs:
- id: bamFile
  type: File
- id: bamIndex
  type: File
- id: referenceFasta
  type: File
- id: referenceFastaFai
  type: File
- id: sample
  type: string
requirements:
- class: InitialWorkDirRequirement
  listing:
  - entryname: example.sh
    entry: 'set -e mkdir -p $(inputs.outputDir) smoove call  --outdir $(inputs.outputDir)  --name
      $(inputs.sample)  --fasta $(inputs.referenceFasta.path)  $(inputs.bamFile.path) '
- class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
- sh
- example.sh
