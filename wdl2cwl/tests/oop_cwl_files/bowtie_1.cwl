class: CommandLineTool
id: Bowtie
inputs:
  - id: readsUpstream
    type:
      - items: File
        type: array
  - id: outputPath
    default: mapped.bam
    type: string
  - id: indexFiles
    type:
      - items: File
        type: array
  - id: best
    default: false
    type: boolean
  - id: strata
    default: false
    type: boolean
  - id: allowContain
    default: false
    type: boolean
  - id: seedmms
    type:
      - int
      - 'null'
  - id: seedlen
    type:
      - int
      - 'null'
  - id: k
    type:
      - int
      - 'null'
  - id: samRG
    type:
      - string
      - 'null'
  - id: picardXmx
    default: 4G
    type: string
  - id: threads
    default: 1
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/mulled-v2-bfe71839265127576d3cd749c056e7b168308d56:1d8bec77b352cdcf3e9ff3d20af238b33ed96eae-0
    type: string
outputs:
  - id: outputBam
    type: File
    outputBinding:
        glob: $(inputs.outputPath)
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-bfe71839265127576d3cd749c056e7b168308d56:1d8bec77b352cdcf3e9ff3d20af238b33ed96eae-0
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
