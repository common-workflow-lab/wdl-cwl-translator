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
outputs:
  - id: smooveVcf
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/$(inputs.sample)-smoove.vcf.gz
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/smoove:0.2.5--0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p $(inputs["outputDir"]===null?"":inputs["outputDir"]["class"]==="File"? "$(inputs.outputDir.path)": "$(inputs.outputDir)")
            smoove call \
            --outdir $(inputs["outputDir"]===null?"":inputs["outputDir"]["class"]==="File"? "$(inputs.outputDir.path)": "$(inputs.outputDir)") \
            --name $(inputs["sample"]===null?"":inputs["sample"]["class"]==="File"? "$(inputs.sample.path)": "$(inputs.sample)") \
            --fasta $(inputs["referenceFasta"]===null?"":inputs["referenceFasta"]["class"]==="File"? "$(inputs.referenceFasta.path)": "$(inputs.referenceFasta)") \
            $(inputs["bamFile"]===null?"":inputs["bamFile"]["class"]==="File"? "$(inputs.bamFile.path)": "$(inputs.bamFile)")
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 15360
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes* 60)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
