class: CommandLineTool
id: Annotate
inputs:
  - id: columns
    default: []
    type:
        items: string
        type: array
  - id: force
    default: false
    type: boolean
  - id: keepSites
    default: false
    type: boolean
  - id: noVersion
    default: false
    type: boolean
  - id: samples
    default: []
    type:
        items: string
        type: array
  - id: singleOverlaps
    default: false
    type: boolean
  - id: removeAnns
    default: []
    type:
        items: string
        type: array
  - id: inputFile
    type: File
  - id: outputPath
    default: output.vcf.gz
    type: string
  - id: annsFile
    type:
      - File
      - 'null'
  - id: collapse
    type:
      - string
      - 'null'
  - id: exclude
    type:
      - string
      - 'null'
  - id: headerLines
    type:
      - File
      - 'null'
  - id: newId
    type:
      - string
      - 'null'
  - id: include
    type:
      - string
      - 'null'
  - id: markSites
    type:
      - string
      - 'null'
  - id: regions
    type:
      - string
      - 'null'
  - id: regionsFile
    type:
      - File
      - 'null'
  - id: renameChrs
    type:
      - File
      - 'null'
  - id: samplesFile
    type:
      - File
      - 'null'
  - id: threads
    default: 0
    type: int
  - id: memory
    default: 256M
    type: string
  - id: dockerImage
    default: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
    type: string
outputs:
  - id: outputVcf
    type: File
    outputBinding:
        glob: $(inputs.outputPath)
  - id: outputVcfIndex
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputPath + ".tbi")
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p "\$(dirname $(inputs.outputPath))"
            bcftools annotate \
            --no-version \
            -o $(inputs.outputPath) \
            -O $(inputs.outputPath.split('/').reverse()[0] !== inputs.outputPath.split('/').reverse()[0].replace(/\.gz$/, '') ? "z" : "v") \
            $(inputs.annsFile === null ? "" : "--annotations " + inputs.annsFile.path) \
            $(inputs.collapse === null ? "" : "--collapse " + inputs.collapse) \
            $(inputs.columns.length > 0 ? "--columns" : "") $(inputs.columns.join(",")) \
            $(inputs.exclude === null ? "" : "--exclude " + inputs.exclude) \
            $(inputs.force ? "--force" : "") \
            $(inputs.headerLines === null ? "" : "--header-lines " + inputs.headerLines.path) \
            $(inputs.newId === null ? "" : "--set-id " + inputs.newId) \
            $(inputs.include === null ? "" : "--include " + inputs.include) \
            $(inputs.keepSites ? "--keep-sites" : "") \
            $(inputs.markSites === null ? "" : "--mark-sites " + inputs.markSites) \
            $(inputs.noVersion ? "--no-version" : "") \
            $(inputs.regions === null ? "" : "--regions " + inputs.regions) \
            $(inputs.regionsFile === null ? "" : "--regions-file " + inputs.regionsFile.path) \
            $(inputs.renameChrs === null ? "" : "--rename-chrs " + inputs.renameChrs.path) \
            $(inputs.samples.length > 0 ? "--samples" : "") $(inputs.samples.join(",")) \
            $(inputs.samplesFile === null ? "" : "--samples-file " + inputs.samplesFile.path) \
            $(inputs.singleOverlaps ? "--single-overlaps" : "") \
            $(inputs.removeAnns.length > 0 ? "--remove" : "") $(inputs.removeAnns.join(",")) \
            $(inputs.inputFile.path)
            # --no-version is for reproducibility
            $(inputs.outputPath.split('/').reverse()[0] !== inputs.outputPath.split('/').reverse()[0].replace(/\.gz$/, '') ? 'bcftools index --tbi ' + inputs.outputPath  : "")
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    ramMin: |-
        ${
        var unit = inputs.memory.match(/[a-zA-Z]+/g).join("");
        var value = parseInt(inputs.memory.match(/[0-9]+/g));
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
