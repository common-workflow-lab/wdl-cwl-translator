cwlVersion: v1.2
$graph:
  - class: CommandLineTool
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
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))"
                bcftools annotate \
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

                $(inputs.outputPath.split('/').reverse()[0] !== inputs.outputPath.split('/').reverse()[0].replace(/\.gz$/, '') ? 'bcftools index --tbi ' + inputs.outputPath : "")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
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
      - script.bash
  - class: CommandLineTool
    id: Sort
    inputs:
      - id: inputFile
        type: File
      - id: outputPath
        default: output.vcf.gz
        type: string
      - id: tmpDir
        default: ./sorting-tmp
        type: string
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
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPath))" $(inputs.tmpDir)
                bcftools sort \
                -o $(inputs.outputPath) \
                -O $(inputs.outputPath.split('/').reverse()[0] !== inputs.outputPath.split('/').reverse()[0].replace(/\.gz$/, '') ? "z" : "v") \
                -T $(inputs.tmpDir) \
                $(inputs.inputFile.path)

                $(inputs.outputPath.split('/').reverse()[0] !== inputs.outputPath.split('/').reverse()[0].replace(/\.gz$/, '') ? 'bcftools index --tbi ' + inputs.outputPath : "")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
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
      - script.bash
  - class: CommandLineTool
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
            glob: '$(inputs.outputPath === null ? inputs.inputVcf.basename + ".stats"
                : inputs.outputPath)'
    requirements:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p \$(dirname $(inputs.outputPath === null ? inputs.inputVcf.basename + ".stats" : inputs.outputPath))
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
                $(inputs.samples.length > 0 ? "--samples" : "") $(inputs.samples.join(",")) \
                $(inputs.samplesFile === null ? "" : "--samples-file " + inputs.samplesFile.path) \
                $(inputs.targets === null ? "" : "--targets " + inputs.targets) \
                $(inputs.targetsFile === null ? "" : "--targets-file " + inputs.targetsFile.path) \
                $(inputs.userTsTv === null ? "" : "--user-tstv " + inputs.userTsTv) \
                --threads $(inputs.threads) \
                $(inputs.verbose ? "--verbose" : "") \
                $(inputs.inputVcf.path) $(inputs.compareVcf === null ? "" : inputs.compareVcf.path) > $(inputs.outputPath === null ? inputs.inputVcf.basename + ".stats" : inputs.outputPath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: $(inputs.threads + 1)
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
      - script.bash
  - class: CommandLineTool
    id: View
    inputs:
      - id: inputFile
        type: File
      - id: outputPath
        default: output.vcf
        type: string
      - id: excludeUncalled
        default: false
        type: boolean
      - id: exclude
        type:
          - string
          - 'null'
      - id: include
        type:
          - string
          - 'null'
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
            glob: $("output.vcf")
      - id: outputVcfIndex
        type:
          - File
          - 'null'
        outputBinding:
            glob: $("output.vcf" + ".tbi")
    requirements:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $("output.vcf"))"
                bcftools view \
                $(inputs.exclude === null ? "" : "--exclude " + inputs.exclude) \
                $(inputs.include === null ? "" : "--include " + inputs.include) \
                $(inputs.excludeUncalled ? "--exclude-uncalled" : "") \
                -o $("output.vcf") \
                -O $(inputs.outputPath.split('/').reverse()[0] !== inputs.outputPath.split('/').reverse()[0].replace(/\.gz$/, '') ? "z" : "v") \
                $(inputs.inputFile.path)

                $(inputs.outputPath.split('/').reverse()[0] !== inputs.outputPath.split('/').reverse()[0].replace(/\.gz$/, '') ? 'bcftools index --tbi ' + "output.vcf" : "")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
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
      - script.bash
