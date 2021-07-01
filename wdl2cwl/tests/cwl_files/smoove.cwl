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
    entry: "\nset -e\nmkdir -p $(inputs.outputDir)\nsmoove call \\\n--outdir $(inputs.outputDir)\
      \ \\\n--name $(inputs.sample) \\\n--fasta $(inputs.referenceFasta.path) \\\n\
      $(inputs.bamFile.path)\n"
- class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
- sh
- example.sh
