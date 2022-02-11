class: CommandLineTool
id: Bowtie
inputs:
  - id: readsUpstream
    type:
        items: File
        type: array
  - id: readsDownstream
    default: []
    type:
        items: File
        type: array
  - id: outputPath
    default: mapped.bam
    type: string
  - id: indexFiles
    type:
        items: File
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
  - id: memory
    type: string
  - id: timeMinutes
    type:
      - int
      - 'null'
  - id: dockerImage
    default: quay.io/biocontainers/mulled-v2-bfe71839265127576d3cd749c056e7b168308d56:1d8bec77b352cdcf3e9ff3d20af238b33ed96eae-0
    type: string
outputs:
  - id: outputBam
    type: File
    outputBinding:
        glob: $(inputs.outputPath)
  - id: outputBamIndex
    type: File
    outputBinding:
        glob: $(inputs.outputPath.replace("\.bam$", ".bai") )
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            set -e -o pipefail
            mkdir -p "\$(dirname $(inputs.outputPath))"
            bowtie \
            -q \
            --sam \
            $(inputs.seedmms === null ? "" : "--seedmms " + inputs.seedmms) \
            $(inputs.seedlen === null ? "" : "--seedlen " + inputs.seedlen) \
            $(inputs.k === null ? "" : "-k " + inputs.k) \
            $(inputs.best ? "--best" : "") \
            $(inputs.strata ? "--strata" : "") \
            $(inputs.allowContain ? "--allow-contain" : "") \
            --threads  $(inputs.threads) \
            $(inputs.samRG === null ? "" : "--sam-RG '" + inputs.samRG)$(inputs.samRG === null ? "" : "'") \
            $(inputs.indexFiles[0].replace("(\.rev)?\.[0-9]\.ebwt$", "") ) \
            $(inputs.readsDownstream.length > 0 ? "-1" : "") $(inputs.readsUpstream.map(function(el) {return el.path}).join(",")) \
            $(inputs.readsDownstream.length > 0 ? "-2" : "") $(inputs.readsDownstream.map(function(el) {return el.path}).join(",")) \
            | picard -Xmx$(inputs.picardXmx) SortSam \
            INPUT=/dev/stdin \
            OUTPUT=$(inputs.outputPath) \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-bfe71839265127576d3cd749c056e7b168308d56:1d8bec77b352cdcf3e9ff3d20af238b33ed96eae-0
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: |-
        ${
        var unit = inputs.memory.match(/[a-zA-Z]+/g).join("");
        var value = parseInt(`${inputs.memory}`.match(/[0-9]+/g));
        var memory = "";
        if(unit==="KiB") memory = value/1024;
        else if(unit==="MiB") memory = value;
        else if(unit==="GiB") memory = value*1024;
        else if(unit==="TiB") memory = value*1024*1024;
        else if(unit==="B") memory = value/(1024*1024);
        else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
        else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
        else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
        else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
        return parseInt(memory);
        }
    outdirMin: 1024
  - class: ToolTimeLimit
    timelimit: $(1 + Math.ceil((function(size_of=0){(function () {var new_array =
        []; [ inputs.readsUpstream, inputs.readsDownstream ].forEach(function(value,
        index, obj) {value.forEach(function(sub_value, sub_index, sub_obj) {new_array.push(sub_value);});});
        return new_array;})().forEach(function(element){ if (element) {size_of +=
        element.size}})}) / 1000^3 * 300 / inputs.threads)  * 60)
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
