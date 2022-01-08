class: CommandLineTool
id: Stats
inputs:
  - id: inputVcf
    type: File
  - id: inputVcfIndex
    type: File
  - id: outputPath
    type:
      - string
      - 'null'
  - id: firstAlleleOnly
    default: false
    type: boolean
  - id: splitByID
    default: false
    type: boolean
  - id: samples
    default: []
    type:
        items: string
        type: array
  - id: verbose
    default: false
    type: boolean
  - id: compareVcf
    type:
      - File
      - 'null'
  - id: compareVcfIndex
    type:
      - File
      - 'null'
  - id: afBins
    type:
      - string
      - 'null'
  - id: afTag
    type:
      - string
      - 'null'
  - id: collapse
    type:
      - string
      - 'null'
  - id: depth
    type:
      - string
      - 'null'
  - id: exclude
    type:
      - string
      - 'null'
  - id: exons
    type:
      - File
      - 'null'
  - id: applyFilters
    type:
      - string
      - 'null'
  - id: fastaRef
    type:
      - File
      - 'null'
  - id: fastaRefIndex
    type:
      - File
      - 'null'
  - id: include
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
  - id: samplesFile
    type:
      - File
      - 'null'
  - id: targets
    type:
      - string
      - 'null'
  - id: targetsFile
    type:
      - File
      - 'null'
  - id: userTsTv
    type:
      - string
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
  - id: stats
    type: File
    outputBinding:
        glob: "$(inputs.outputPath == null ? inputs.inputVcf.basename + '.stats' :\
            \ inputs.outputPath)"
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p \$(dirname $(inputs.outputPath == null ? inputs.inputVcf.basename + '.stats' : inputs.outputPath))
            bcftools stats \
            $(inputs.afBins === null ? "" : "--af-bins " + inputs.afBins) \
            $(inputs.afTag === null ? "" : "--af-tag " + inputs.afTag) \
            $(inputs.firstAlleleOnly ? "--1st-allele-only" : "") \
            $(inputs.collapse === null ? "" : "--collapse " + inputs.collapse) \
            $(inputs.depth === null ? "" : "--depth " + inputs.depth) \
            $(inputs.exclude === null ? "" : "--exclude " + inputs.exclude) \
            $(inputs.exons === null ? "" : "--exons " + inputs.exons.path) \
            $(inputs.applyFilters === null ? "" : "--apply-filters " + inputs.applyFilters) \
            $(inputs.fastaRef === null ? "" : "--fasta-ref " + inputs.fastaRef.path) \
            $(inputs.include === null ? "" : "--include " + inputs.include) \
            $(inputs.splitByID ? "--split-by-ID" : "") \
            $(inputs.regions === null ? "" : "--regions " + inputs.regions) \
            $(inputs.regionsFile === null ? "" : "--regions-file " + inputs.regionsFile.path) \
            $(inputs.samples.length > 0 ? "--samples" : "") $(inputs.samples.map(function(el) {return el.path}).join(",")) \
            $(inputs.samplesFile === null ? "" : "--samples-file " + inputs.samplesFile.path) \
            $(inputs.targets === null ? "" : "--targets " + inputs.targets) \
            $(inputs.targetsFile === null ? "" : "--targets-file " + inputs.targetsFile.path) \
            $(inputs.userTsTv === null ? "" : "--user-tstv " + inputs.userTsTv) \
            --threads $(inputs.threads) \
            $(inputs.verbose ? "--verbose" : "") \
            $(inputs.inputVcf.path) $(inputs.compareVcf == null ? "" : inputs.compareVcf.path) > $(inputs.outputPath == null ? inputs.inputVcf.basename + '.stats' : inputs.outputPath)
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    coresMin: $(inputs.threads + 1)
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
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
