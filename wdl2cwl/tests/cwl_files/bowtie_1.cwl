class: CommandLineTool
id: Bowtie
inputs:
  - id: readsUpstream
    type:
      - items: File
        type: array
  - id: indexFiles
    type:
      - items: File
        type: array
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
  - id: outputPath
    default: mapped.bam
    type: string
  - id: best
    default: false
    type: boolean
  - id: strata
    default: false
    type: boolean
  - id: allowContain
    default: false
    type: boolean
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
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e -o pipefail
            mkdir -p "\$(dirname $(inputs.outputPath))"
            bowtie \
            -q \
            --sam \
            $(inputs.seedmms === null ? "" : "--seedmms " + inputs.seedmms ) \
            $(inputs.seedlen === null ? "" : "--seedlen " + inputs.seedlen ) \
            $(inputs.k === null ? "" : "-k " + inputs.k ) \
            $(inputs.best ? "--best" : "") \
            $(inputs.strata ? "--strata" : "") \
            $(inputs.allowContain ? "--allow-contain" : "") \
            --threads $(inputs.threads) \
            $(inputs.samRG === null ? "" : "--sam-RG '" + inputs.samRG )$(inputs.samRG === null  ? "" : "'") \
            $(inputs.indexFiles[0].replace("(\.rev)?\.[0-9]\.ebwt$","")) \
            $(inputs.readsUpstream.map(function(el) { return el.path}).join(",")) \
            | picard -Xmx$(inputs.picardXmx) SortSam \
            INPUT=/dev/stdin \
            OUTPUT=$(inputs.outputPath) \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
