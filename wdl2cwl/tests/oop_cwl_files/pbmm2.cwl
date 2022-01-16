class: CommandLineTool
id: Mapping
inputs:
  - id: presetOption
    type: string
  - id: sort
    default: true
    type: boolean
  - id: sample
    type: string
  - id: referenceMMI
    type: File
  - id: queryFile
    type: File
  - id: cores
    default: 4
    type: int
  - id: memory
    default: 30G
    type: string
  - id: dockerImage
    default: quay.io/biocontainers/pbmm2:1.3.0--h56fc30b_1
    type: string
outputs:
  - id: outputAlignmentFile
    type: File
    outputBinding:
        glob: $(inputs.sample + ".align.bam")
  - id: outputIndexFile
    type: File
    outputBinding:
        glob: $(inputs.sample + ".align.bam.bai")
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/pbmm2:1.3.0--h56fc30b_1
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            pbmm2 align \
            --preset $(inputs.presetOption) \
            $(inputs.sort ? "--sort" : "") \
            -j $(inputs.cores) \
            $(inputs.referenceMMI.path) \
            $(inputs.queryFile.path) \
            --sample $(inputs.sample) \
            $(inputs.sample).align.bam
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    coresMin: $(inputs.cores)
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
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
