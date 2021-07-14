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
outputs:
  - id: smooveVcf
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/$(inputs.sample)-smoove.vcf.gz
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
