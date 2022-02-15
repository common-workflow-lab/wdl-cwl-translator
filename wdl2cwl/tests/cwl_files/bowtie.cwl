cwlVersion: v1.2
id: Bowtie
class: CommandLineTool
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
            $("--threads " + inputs.threads) \
            $(inputs.samRG === null ? "" : "--sam-RG '" + inputs.samRG)$(inputs.samRG === null ? "" : "'") \
            $(inputs.indexFiles[0].replace("(\\.rev)?\\.[0-9]\\.ebwt$", "") ) \
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
        else throw "Unknown units: " + unit;
        return parseInt(memory);
        }
    outdirMin: 1024
  - class: ToolTimeLimit
    timelimit: $(1 + Math.ceil((function(size_of=0){(function () {var new_array =
        []; [ inputs.readsUpstream, inputs.readsDownstream ].forEach(function(value,
        index, obj) {value.forEach(function(sub_value, sub_index, sub_obj) {new_array.push(sub_value);});});
        return new_array;})().forEach(function(element){ if (element) {size_of +=
        element.size}})}) / 1000^3 * 300 / inputs.threads)  * 60)
inputs:
  - id: readsUpstream
    doc: The first-/single-end fastq files.
    type:
        items: File
        type: array
  - id: readsDownstream
    doc: The second-end fastq files.
    default: []
    type:
        items: File
        type: array
  - id: outputPath
    doc: The location the output BAM file should be written to.
    default: mapped.bam
    type: string
  - id: indexFiles
    doc: The index files for bowtie.
    type:
        items: File
        type: array
  - id: best
    doc: Equivalent to bowtie's `--best` flag.
    default: false
    type: boolean
  - id: strata
    doc: Equivalent to bowtie's `--strata` flag.
    default: false
    type: boolean
  - id: allowContain
    doc: Equivalent to bowtie's `--allow-contain` flag.
    default: false
    type: boolean
  - id: seedmms
    doc: Equivalent to bowtie's `--seedmms` option.
    type:
      - int
      - 'null'
  - id: seedlen
    doc: Equivalent to bowtie's `--seedlen` option.
    type:
      - int
      - 'null'
  - id: k
    doc: Equivalent to bowtie's `-k` option.
    type:
      - int
      - 'null'
  - id: samRG
    doc: Equivalent to bowtie's `--sam-RG` option.
    type:
      - string
      - 'null'
  - id: picardXmx
    doc: The maximum memory available to the picard (used for sorting the output).
        Should be lower than `memory` to accommodate JVM overhead and bowtie's memory
        usage.
    default: 4G
    type: string
  - id: threads
    doc: The number of threads to use.
    default: 1
    type: int
  - id: memory
    doc: The amount of memory this job will use.
    type: string
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    type:
      - int
      - 'null'
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/mulled-v2-bfe71839265127576d3cd749c056e7b168308d56:1d8bec77b352cdcf3e9ff3d20af238b33ed96eae-0
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: outputBam
    doc: Output alignment file.
    type: File
    outputBinding:
        glob: $(inputs.outputPath)
  - id: outputBamIndex
    doc: Index of output alignment file.
    type: File
    outputBinding:
        glob: $(inputs.outputPath.replace("\\.bam$", ".bai") )
