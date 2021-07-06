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
  - id: outputDir
    default: ./smoove
    type: string
  - id: memory
    default: 15G
    type: string
  - id: timeMinutes
    default: 1440
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/smoove:0.2.5--0
    type: string
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p $(inputs.outputDir)
            smoove call \
            --outdir $(inputs.outputDir) \
            --name $(inputs.sample) \
            --fasta $(inputs.referenceFasta.path) \
            $(inputs.bamFile.path)
  - class: InlineJavascriptRequirement
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
