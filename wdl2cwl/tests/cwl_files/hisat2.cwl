class: CommandLineTool
id: Hisat2
inputs:
  - id: inputR1
    type: File
  - id: inputR2
    type:
      - File
      - 'null'
  - id: indexFiles
    type:
        items: File
        type: array
  - id: outputBam
    type: string
  - id: sample
    type: string
  - id: library
    type: string
  - id: readgroup
    type: string
  - id: platform
    default: illumina
    type: string
  - id: downstreamTranscriptomeAssembly
    default: true
    type: boolean
  - id: summaryFilePath
    type:
      - string
      - 'null'
  - id: sortMemoryPerThreadGb
    default: 2
    type: int
  - id: compressionLevel
    default: 1
    type: int
  - id: sortThreads
    type:
      - int
      - 'null'
  - id: threads
    default: 4
    type: int
  - id: memoryGb
    type:
      - int
      - 'null'
  - id: timeMinutes
    type:
      - int
      - 'null'
  - id: dockerImage
    default: quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0
    type: string
outputs:
  - id: bamFile
    type: File
    outputBinding:
        glob: $(inputs.outputBam)
  - id: summaryFile
    type: File
    outputBinding:
        glob: "$(inputs.summaryFilePath === null ? inputs.outputBam.split('/').reverse()[0].replace(/\\\
            .bam$/, '') + \".summary.txt\" : inputs.summaryFilePath)"
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            set -e -o pipefail
            mkdir -p "\$(dirname $(inputs.outputBam))"
            hisat2 \
            -p $(inputs.threads) \
            -x $(inputs.indexFiles[0].replace("\.[0-9]\.ht2", "") ) \
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
            -@  $(inputs.totalSortThreads) \
            -m $(inputs.sortMemoryPerThreadGb)G \
            -l $(inputs.compressionLevel) \
            - \
            -o $(inputs.outputBam)
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: |-
        ${
        var unit = "G";
        var value = parseInt(`${[inputs.memoryGb, 1 + Math.ceil((function(size_of=0){inputs.indexFiles.forEach(function(element){ if (element) {size_of += element.size}})}) / 1024^3*1.2)  + inputs.sortMemoryPerThreadGb*[inputs.sortThreads, inputs.threads === 1 ? 1 : 1 + Math.ceil(inputs.threads/4.0) ].find(element => element !== null) ].find(element => element !== null) }`.match(/[0-9]+/g));
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
    timelimit: '$(1 + Math.ceil((function(size_of=0){[inputs.inputR1.path, inputs.inputR2
        === null ? "" : inputs.inputR2.path].forEach(function(element){ if (element)
        {size_of += element.size}})}) / 1024^3*180/inputs.threads)  * 60)'
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
