cwlVersion: v1.2
$graph:
  - class: CommandLineTool
    id: Annotate
    inputs:
      - id: columns
        doc: Comma-separated list of columns or tags to carry over from the annotation
            file (see man page for details).
        default: []
        type:
            items: string
            type: array
      - id: force
        doc: Continue even when parsing errors, such as undefined tags, are encountered.
        default: false
        type: boolean
      - id: keepSites
        doc: Keep sites which do not pass -i and -e expressions instead of discarding
            them.
        default: false
        type: boolean
      - id: noVersion
        doc: Do not append version and command line information to the output VCF
            header.
        default: false
        type: boolean
      - id: samples
        doc: List of samples for sample stats, "-" to include all samples.
        default: []
        type:
            items: string
            type: array
      - id: singleOverlaps
        doc: keep memory requirements low with very large annotation files.
        default: false
        type: boolean
      - id: removeAnns
        doc: List of annotations to remove (see man page for details).
        default: []
        type:
            items: string
            type: array
      - id: inputFile
        doc: A vcf or bcf file.
        type: File
      - id: inputFileIndex
        doc: The index for the input vcf or bcf.
        type:
          - File
          - 'null'
      - id: outputPath
        doc: The location the output VCF file should be written.
        default: output.vcf.gz
        type: string
      - id: annsFile
        doc: Bgzip-compressed and tabix-indexed file with annotations (see man page
            for details).
        type:
          - File
          - 'null'
      - id: annsFileIndex
        doc: The index for annsFile.
        type:
          - File
          - 'null'
      - id: collapse
        doc: Treat as identical records with <snps|indels|both|all|some|none>, see
            man page for details.
        type:
          - string
          - 'null'
      - id: exclude
        doc: Exclude sites for which the expression is true (see man page for details).
        type:
          - string
          - 'null'
      - id: headerLines
        doc: Lines to append to the VCF header (see man page for details).
        type:
          - File
          - 'null'
      - id: newId
        doc: Assign ID on the fly (e.g. --set-id +'%CHROM\_%POS').
        type:
          - string
          - 'null'
      - id: include
        doc: Select sites for which the expression is true (see man page for details).
        type:
          - string
          - 'null'
      - id: markSites
        doc: Annotate sites which are present ('+') or absent ('-') in the -a file
            with a new INFO/TAG flag.
        type:
          - string
          - 'null'
      - id: regions
        doc: Restrict to comma-separated list of regions.
        type:
          - string
          - 'null'
      - id: regionsFile
        doc: Restrict to regions listed in a file.
        type:
          - File
          - 'null'
      - id: renameChrs
        doc: rename chromosomes according to the map in file (see man page for details).
        type:
          - File
          - 'null'
      - id: samplesFile
        doc: File of samples to include.
        type:
          - File
          - 'null'
      - id: threads
        doc: Number of extra decompression threads [0].
        default: 0
        type: int
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
        type: string
    outputs:
      - id: outputVcf
        doc: Annotated VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVcfIndex
        doc: Index of the annotated VCF file.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
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
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
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
      - class: ToolTimeLimit
        timelimit: $(60 + Math.ceil((function(size_of=0){inputs.inputFile.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: Filter
    inputs:
      - id: vcf
        doc: The VCF file to operate on.
        type: File
      - id: vcfIndex
        doc: The index for the VCF file.
        type: File
      - id: include
        doc: Equivalent to the `-i` option.
        type:
          - string
          - 'null'
      - id: exclude
        type:
          - string
          - 'null'
      - id: softFilter
        type:
          - string
          - 'null'
      - id: outputPath
        doc: The location the output VCF file should be written.
        default: ./filtered.vcf.gz
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 256M
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
        type: string
    outputs:
      - id: outputVcf
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e 
                mkdir -p "\$(dirname $(inputs.outputPath))"
                bcftools \
                filter \
                $(inputs.include === null ? "" : "-i " + inputs.include) \
                $(inputs.exclude === null ? "" : "-e " + inputs.exclude) \
                $(inputs.softFilter === null ? "" : "-s " + inputs.softFilter) \
                $(inputs.vcf.path) \
                -O z \
                -o $(inputs.outputPath)
                bcftools index --tbi $(inputs.outputPath)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
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
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.vcf.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: Sort
    inputs:
      - id: inputFile
        doc: A vcf or bcf file.
        type: File
      - id: outputPath
        doc: The location the output VCF file should be written.
        default: output.vcf.gz
        type: string
      - id: tmpDir
        doc: The location of the temporary files during the bcftools sorting.
        default: ./sorting-tmp
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 256M
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
        type: string
    outputs:
      - id: outputVcf
        doc: Sorted VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
      - id: outputVcfIndex
        doc: Index of sorted VCF file.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputPath + ".tbi")
    requirements:
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
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
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
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.inputFile.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: Stats
    inputs:
      - id: inputVcf
        doc: The VCF to be analysed.
        type: File
      - id: inputVcfIndex
        doc: The index for the input VCF.
        type: File
      - id: outputPath
        doc: The location the output VCF file should be written.
        type:
          - string
          - 'null'
      - id: firstAlleleOnly
        doc: Include only 1st allele at multiallelic sites.
        default: false
        type: boolean
      - id: splitByID
        doc: Collect stats for sites with ID separately (known vs novel).
        default: false
        type: boolean
      - id: samples
        doc: List of samples for sample stats, "-" to include all samples.
        default: []
        type:
            items: string
            type: array
      - id: verbose
        doc: Produce verbose per-site and per-sample output.
        default: false
        type: boolean
      - id: compareVcf
        doc: When inputVcf and compareVCF are given, the program generates separate
            stats for intersection and the complements. By default only sites are
            compared, samples must be given to include also sample columns.
        type:
          - File
          - 'null'
      - id: compareVcfIndex
        doc: Index for the compareVcf.
        type:
          - File
          - 'null'
      - id: afBins
        doc: |-
            Allele frequency bins, a list (0.1,0.5,1) or a file (0.1
            0.5
            1).
        type:
          - string
          - 'null'
      - id: afTag
        doc: Allele frequency tag to use, by default estimated from AN,AC or GT.
        type:
          - string
          - 'null'
      - id: collapse
        doc: Treat as identical records with <snps|indels|both|all|some|none>, see
            man page for details.
        type:
          - string
          - 'null'
      - id: depth
        doc: 'Depth distribution: min,max,bin size [0,500,1].'
        type:
          - string
          - 'null'
      - id: exclude
        doc: Exclude sites for which the expression is true (see man page for details).
        type:
          - string
          - 'null'
      - id: exons
        doc: Tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based,
            inclusive, bgzip compressed).
        type:
          - File
          - 'null'
      - id: applyFilters
        doc: Require at least one of the listed FILTER strings (e.g. "PASS,.").
        type:
          - string
          - 'null'
      - id: fastaRef
        doc: Faidx indexed reference sequence file to determine INDEL context.
        type:
          - File
          - 'null'
      - id: fastaRefIndex
        doc: Index file (.fai) for fastaRef. Must be supplied if fastaRef is supplied.
        type:
          - File
          - 'null'
      - id: include
        doc: Select sites for which the expression is true (see man page for details).
        type:
          - string
          - 'null'
      - id: regions
        doc: Restrict to comma-separated list of regions.
        type:
          - string
          - 'null'
      - id: regionsFile
        doc: Restrict to regions listed in a file.
        type:
          - File
          - 'null'
      - id: samplesFile
        doc: File of samples to include.
        type:
          - File
          - 'null'
      - id: targets
        doc: Similar to regions but streams rather than index-jumps.
        type:
          - string
          - 'null'
      - id: targetsFile
        doc: Similar to regionsFile but streams rather than index-jumps.
        type:
          - File
          - 'null'
      - id: userTsTv
        doc: <TAG[:min:max:n]>. Collect Ts/Tv stats for any tag using the given binning
            [0:1:100].
        type:
          - string
          - 'null'
      - id: threads
        doc: Number of extra decompression threads [0].
        default: 0
        type: int
      - id: memory
        doc: The amount of memory this job will use.
        default: 256M
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
        type: string
    outputs:
      - id: stats
        doc: Text file stats which is suitable for machine processing and can be plotted
            using plot-vcfstats.
        type: File
        outputBinding:
            glob: '$(inputs.outputPath === null ? inputs.inputVcf.basename + ".stats"
                : inputs.outputPath)'
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p \$(dirname $(inputs.outputPath === null ? inputs.inputVcf.basename + ".stats" : inputs.outputPath))
                mkdir fastaRef_dir  # to ensure correct localization
                ln -s $(inputs.fastaRef === null ? "" : inputs.fastaRef.path) fastaRef_dir/\$(basename $(inputs.fastaRef === null ? "" : inputs.fastaRef.path))
                ln -s $(inputs.fastaRefIndex === null ? "" : inputs.fastaRefIndex.path) fastaRef_dir/\$(basename $(inputs.fastaRefIndex === null ? "" : inputs.fastaRefIndex.path))
                bcftools stats \
                $(inputs.afBins === null ? "" : "--af-bins " + inputs.afBins) \
                $(inputs.afTag === null ? "" : "--af-tag " + inputs.afTag) \
                $(inputs.firstAlleleOnly ? "--1st-allele-only" : "") \
                $(inputs.collapse === null ? "" : "--collapse " + inputs.collapse) \
                $(inputs.depth === null ? "" : "--depth " + inputs.depth) \
                $(inputs.exclude === null ? "" : "--exclude " + inputs.exclude) \
                $(inputs.exons === null ? "" : "--exons " + inputs.exons.path) \
                $(inputs.applyFilters === null ? "" : "--apply-filters " + inputs.applyFilters) \
                $(inputs.fastaRef === null ? "" : "--fasta-ref fastaRef_dir/$(basename " + inputs.fastaRef.path + ")") \
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
                sed -i "s=\$(dirname $(inputs.inputVcf.path))/==g" $(inputs.outputPath === null ? inputs.inputVcf.basename + ".stats" : inputs.outputPath)  # for reproducibility
                sed -i "s=\$(dirname $(inputs.fastaRef === null ? "" : inputs.fastaRef.path))/==g" $(inputs.outputPath === null ? inputs.inputVcf.basename + ".stats" : inputs.outputPath)  # for reproducibility
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
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
      - class: ToolTimeLimit
        timelimit: '$(1 + 2 * Math.ceil((function(size_of=0){[inputs.inputVcf.path,
            inputs.compareVcf === null ? "" : inputs.compareVcf.path].filter(function(element)
            { return element !== null }) .forEach(function(element){ if (element)
            {size_of += element.size}})}) / 1000^3)  * 60)'
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: View
    inputs:
      - id: inputFile
        doc: A vcf or bcf file.
        type: File
      - id: outputPath
        doc: The location the output VCF file should be written.
        default: output.vcf
        type: string
      - id: excludeUncalled
        doc: Exclude sites without a called genotype (see man page for details).
        default: false
        type: boolean
      - id: exclude
        doc: Exclude sites for which the expression is true (see man page for details).
        type:
          - string
          - 'null'
      - id: include
        doc: Select sites for which the expression is true (see man page for details).
        type:
          - string
          - 'null'
      - id: memory
        doc: The amount of memory this job will use.
        default: 256M
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
        type: string
    outputs:
      - id: outputVcf
        doc: VCF file.
        type: File
        outputBinding:
            glob: $("output.vcf")
      - id: outputVcfIndex
        doc: Index of VCF file.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $("output.vcf" + ".tbi")
    requirements:
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
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2
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
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.inputFile.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
