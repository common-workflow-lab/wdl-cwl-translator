cwlVersion: v1.2
id: Hisat2
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            set -e -o pipefail
            mkdir -p "\$(dirname $(inputs.outputBam))"
            hisat2 \
            -p $(inputs.threads) \
            -x $(inputs.indexFiles[0].replace("\\.[0-9]\\.ht2", "") ) \
            $(inputs.inputR2 === null ? "-U" : "-1") $(inputs.inputR1.path) \
            $(inputs.inputR2 === null ? "" : "-2" + inputs.inputR2.path) \
            --rg-id $(inputs.readgroup) \
            --rg 'SM:$(inputs.sample)' \
            --rg 'LB:$(inputs.library)' \
            --rg 'PL:$(inputs.platform)' \
            $(inputs.downstreamTranscriptomeAssembly ? "--dta" : "") \
            --new-summary \
            --summary-file $(inputs.summaryFilePath === null ? inputs.outputBam.split('/').reverse()[0].replace(/\.bam$/, '') + ".summary.txt" : inputs.summaryFilePath) \
            | samtools sort \
            $("-@ " + inputs.totalSortThreads) \
            -m $(inputs.sortMemoryPerThreadGb)G \
            -l $(inputs.compressionLevel) \
            - \
            -o $(inputs.outputBam)
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: |-
        ${
        var unit = "G";
        var value = parseInt(`${[inputs.memoryGb, 1 + Math.ceil((function(size_of=0){inputs.indexFiles.forEach(function(element){ if (element) {size_of += element.size}})}) / 1000^3 * 1.2)  + inputs.sortMemoryPerThreadGb * [inputs.sortThreads, inputs.threads === 1 ? 1 : 1 + Math.ceil(inputs.threads / 4.0) ].find(function(element) { return element !== null }) ].find(function(element) { return element !== null }) }`.match(/[0-9]+/g));
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
    timelimit: '$(1 + Math.ceil((function(size_of=0){[inputs.inputR1.path, inputs.inputR2
        === null ? "" : inputs.inputR2.path].forEach(function(element){ if (element)
        {size_of += element.size}})}) / 1000^3 * 180 / inputs.threads)  * 60)'
inputs:
  - id: inputR1
    doc: The first-/single-end FastQ file.
    type: File
  - id: inputR2
    doc: The second-end FastQ file.
    type:
      - File
      - 'null'
  - id: indexFiles
    doc: The hisat2 index files.
    type:
        items: File
        type: array
  - id: outputBam
    doc: The location the output BAM file should be written to.
    type: string
  - id: sample
    doc: The sample id.
    type: string
  - id: library
    doc: The library id.
    type: string
  - id: readgroup
    doc: The readgroup id.
    type: string
  - id: platform
    doc: The platform used for sequencing.
    default: illumina
    type: string
  - id: downstreamTranscriptomeAssembly
    doc: Equivalent to hisat2's `--dta` flag.
    default: true
    type: boolean
  - id: summaryFilePath
    doc: Where the summary file should be written.
    type:
      - string
      - 'null'
  - id: sortMemoryPerThreadGb
    doc: The amount of memory for each sorting thread in gigabytes.
    default: 2
    type: int
  - id: compressionLevel
    doc: The compression level of the output BAM.
    default: 1
    type: int
  - id: sortThreads
    doc: The number of threads to use for sorting.
    type:
      - int
      - 'null'
  - id: threads
    doc: The number of threads to use.
    default: 4
    type: int
  - id: memoryGb
    doc: The amount of memory this job will use in gigabytes.
    type:
      - int
      - 'null'
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    type:
      - int
      - 'null'
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: bamFile
    doc: Output BAM file.
    type: File
    outputBinding:
        glob: $(inputs.outputBam)
  - id: summaryFile
    doc: Alignment summary file.
    type: File
    outputBinding:
        glob: "$(inputs.summaryFilePath === null ? inputs.outputBam.split('/').reverse()[0].replace(/\\\
            .bam$/, '') + \".summary.txt\" : inputs.summaryFilePath)"
