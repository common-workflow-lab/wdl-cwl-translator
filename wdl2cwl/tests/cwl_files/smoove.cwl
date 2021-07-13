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
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/smoove:0.2.5--0
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
  - class: ResourceRequirement
    ramMin: 15360
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes* 60)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
