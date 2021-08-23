class: CommandLineTool
id: Hisat2
inputs:
  - id: inputR1
    type: File
  - id: inputR2
    default: ''
    type:
      - File
      - 'null'
  - id: indexFiles
    type:
      - items: File
        type: array
  - id: outputBam
    type: string
  - id: sample
    type: string
  - id: library
    type: string
  - id: readgroup
    type: string
  - id: sortThreads
    default: ''
    type:
      - int
      - 'null'
  - id: memoryGb
    default: ''
    type:
      - int
      - 'null'
  - id: platform
    default: illumina
    type: string
  - id: downstreamTranscriptomeAssembly
    default: true
    type: boolean
  - id: sortMemoryPerThreadGb
    default: 2
    type: int
  - id: compressionLevel
    default: 1
    type: int
  - id: threads
    default: 4
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0
    type: string
outputs:
  - id: bamFile
    type: File
    outputBinding:
        glob: $(inputs.outputBam)
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e -o pipefail
            mkdir -p "\$(dirname $(inputs.outputBam))"
            hisat2 \
            -p $(inputs.threads) \
            -x $(inputs.indexFiles[0]) \
            $(inputs["inputR2"] === "" ? "-U": "-1") $(inputs.inputR1.path) \
            -2$(inputs.inputR2) \
            --rg-id $(inputs.readgroup) \
            --rg 'SM:$(inputs.sample)' \
            --rg 'LB:$(inputs.library)' \
            --rg 'PL:$(inputs.platform)' \
            $(inputs["downstreamTranscriptomeAssembly"] ? "--dta" : "") \
            --new-summary \
            --summary-file $(inputs.summaryFilePath) \
            | samtools sort \
            "-@ " \
            -m $(inputs.sortMemoryPerThreadGb)G \
            -l $(inputs.compressionLevel) \
            - \
            -o $(inputs.outputBam)
  - class: InlineJavascriptRequirement
  - class: ToolTimeLimit
    timelimit: 60
  - class: ResourceRequirement
    coresMin: 2
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
