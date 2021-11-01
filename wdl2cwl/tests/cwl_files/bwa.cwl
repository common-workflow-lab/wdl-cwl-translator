class: CommandLineTool
id: Mem
inputs:
  - id: read1
    type: File
  - id: read2
    type:
      - File
      - 'null'
  - id: outputPrefix
    type: string
  - id: readgroup
    type:
      - string
      - 'null'
  - id: sortThreads
    type:
      - int
      - 'null'
  - id: memoryGb
    type:
      - int
      - 'null'
  - id: sixtyFour
    default: false
    type: boolean
  - id: usePostalt
    default: false
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
    default: quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0
    type: string
outputs:
  - id: outputBam
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix).aln.bam
  - id: outputHla
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputPrefix).hla
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p "\$(dirname $(inputs.outputPrefix))"
            bwa mem \
              -t $(inputs.threads) \
              $(inputs.readgroup === null ? "" : "-R '" + inputs.readgroup )$(inputs.readgroup === null  ? "" : "'") \
              $(inputs.bwaIndex.fastaFile) \
              $(inputs.read1.path) \
              $(inputs.read2 === null ? "" : inputs.read2.path) \
              2> $(inputs.outputPrefix).log.bwamem | \
              $(inputs.usePostalt ? "" : "#") bwa-postalt.js -p $(inputs.outputPrefix).hla $(inputs.bwaIndex.fastaFile)$(inputs.sixtyFour ? ".64.alt" : ".alt") | \
              samtools sort \
              -m $(inputs.sortMemoryPerThreadGb)G \
              -l $(inputs.compressionLevel) \
              - \
              -o $(inputs.outputPrefix).aln.bam
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    ramMin: ${[inputs.memoryGb].find(mem => mem !== null)}G
  - class: ToolTimeLimit
    timelimit: ''
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
