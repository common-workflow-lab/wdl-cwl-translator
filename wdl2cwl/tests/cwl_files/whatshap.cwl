cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: Phase
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                whatshap phase \
                $(inputs.vcf.path) \
                $(inputs.phaseInput.path) \
                $(inputs.outputVCF ? "--output " + "\"" + inputs.outputVCF + "\"" : "") \
                $(inputs.reference ? inputs.reference === null ? "" : "--reference " + "\"" + inputs.reference.path + "\"" : "") \
                $(inputs.tag ? inputs.tag === null ? "" : "--tag " + "\"" + inputs.tag + "\"" : "") \
                $(inputs.algorithm ? inputs.algorithm === null ? "" : "--algorithm " + "\"" + inputs.algorithm + "\"" : "") \
                $(inputs.indels === null ? "" : "--indels") \
                $(inputs.sample ? inputs.sample === null ? "" : "--sample " + "\"" + inputs.sample + "\"" : "") \
                $(inputs.chromosome ? inputs.chromosome === null ? "" : "--chromosome " + "\"" + inputs.chromosome + "\"" : "") \
                $(inputs.threshold ? inputs.threshold === null ? "" : "--threshold " + "\"" + inputs.threshold + "\"" : "") \
                $(inputs.ped ? inputs.ped === null ? "" : "--ped " + "\"" + inputs.ped + "\"" : "")

                tabix -p vcf $(inputs.outputVCF)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: outputVCF
        doc: Output VCF file. Add .gz to the file name to get compressed output. If
            omitted, use standard output.
        type: string
      - id: vcf
        doc: VCF or BCF file with variants to be phased (can be gzip-compressed).
        type: File
      - id: vcfIndex
        doc: Index for the VCF or BCF file with variants to be phased.
        type: File
      - id: phaseInput
        doc: BAM, CRAM, VCF or BCF file(s) with phase information, either through
            sequencing reads (BAM, CRAM) or through phased blocks (VCF, BCF).
        type: File
      - id: phaseInputIndex
        doc: Index of BAM, CRAM, VCF or BCF file(s) with phase information.
        type: File
      - id: reference
        doc: Reference file. Provide this to detect alleles through re-alignment.
            If no index (.fai) exists, it will be created.
        type:
          - File
          - 'null'
      - id: referenceIndex
        doc: Index of reference file.
        type:
          - File
          - 'null'
      - id: tag
        doc: 'Store phasing information with PS tag (standardized) or HP tag (used
            by GATK ReadBackedPhasing) (default: {description: PS).'
        type:
          - string
          - 'null'
      - id: algorithm
        doc: 'Phasing algorithm to use (default: {description: whatshap).'
        type:
          - string
          - 'null'
      - id: indels
        doc: 'Also phase indels (default: {description: do not phase indels).'
        type:
          - boolean
          - 'null'
      - id: sample
        doc: Name of a sample to phase. If not given, all samples in the input VCF
            are phased. Can be used multiple times.
        type:
          - string
          - 'null'
      - id: chromosome
        doc: Name of chromosome to phase. If not given, all chromosomes in the input
            VCF are phased. Can be used multiple times.
        type:
          - string
          - 'null'
      - id: threshold
        doc: 'The threshold of the ratio between the probabilities that a pair of
            reads come from the same haplotype and different haplotypes in the read
            merging model (default: {description: 1000000).'
        type:
          - string
          - 'null'
      - id: ped
        doc: Use pedigree information in PED file to improve phasing (switches to
            PedMEC algorithm). Columns 2, 3, 4 must refer to child, mother, and father
            sample names as used in the VCF and BAM/CRAM. Other columns are ignored.
        type:
          - string
          - 'null'
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: phasedVCF
        doc: VCF file containing phased variants.
        type: File
        outputBinding:
            glob: $(inputs.outputVCF)
      - id: phasedVCFIndex
        doc: Index of phased VCF file.
        type: File
        outputBinding:
            glob: $(inputs.outputVCF + ".tbi")
  - cwlVersion: v1.2
    id: Stats
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                whatshap stats \
                $(inputs.vcf.path) \
                $(inputs.gtf ? inputs.gtf === null ? "" : "--gtf " + "\"" + inputs.gtf + "\"" : "") \
                $(inputs.sample ? inputs.sample === null ? "" : "--sample " + "\"" + inputs.sample + "\"" : "") \
                $(inputs.tsv ? inputs.tsv === null ? "" : "--tsv " + "\"" + inputs.tsv + "\"" : "") \
                $(inputs.blockList ? inputs.blockList === null ? "" : "--block-list " + "\"" + inputs.blockList + "\"" : "") \
                $(inputs.chromosome ? inputs.chromosome === null ? "" : "--chromosome " + "\"" + inputs.chromosome + "\"" : "")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: vcf
        doc: Phased VCF file.
        type: File
      - id: gtf
        doc: Write phased blocks to GTF file.
        type:
          - string
          - 'null'
      - id: sample
        doc: Name of the sample to process. If not given, use first sample found in
            VCF.
        type:
          - string
          - 'null'
      - id: tsv
        doc: Filename to write statistics to (tab-separated).
        type:
          - string
          - 'null'
      - id: blockList
        doc: Filename to write list of all blocks to (one block per line).
        type:
          - string
          - 'null'
      - id: chromosome
        doc: Name of chromosome to process. If not given, all chromosomes in the input
            VCF are considered.
        type:
          - string
          - 'null'
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: phasedGTF
        doc: Phasing statistics for a single VCF file.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.gtf)
      - id: phasedTSV
        doc: Statistics in a tab-separated value format.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.tsv)
      - id: phasedBlockList
        doc: List of the total number of phase sets/blocks.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.blockList)
  - cwlVersion: v1.2
    id: Haplotag
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                whatshap haplotag \
                $(inputs.vcf.path) \
                $(inputs.alignments.path) \
                $(inputs.outputFile ? "--output " + "\"" + inputs.outputFile + "\"" : "") \
                $(inputs.reference ? inputs.reference === null ? "" : "--reference " + "\"" + inputs.reference.path + "\"" : "") \
                $(inputs.regions ? inputs.regions === null ? "" : "--regions " + "\"" + inputs.regions + "\"" : "") \
                $(inputs.sample ? inputs.sample === null ? "" : "--sample " + "\"" + inputs.sample + "\"" : "")

                python3 -c "import pysam; pysam.index('$(inputs.outputFile)')"
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0
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
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: vcf
        doc: VCF file with phased variants (must be gzip-compressed and indexed).
        type: File
      - id: vcfIndex
        doc: Index for the VCF or BCF file with variants to be phased.
        type: File
      - id: alignments
        doc: File (BAM/CRAM) with read alignments to be tagged by haplotype.
        type: File
      - id: alignmentsIndex
        doc: Index for the alignment file.
        type: File
      - id: outputFile
        doc: Output file. If omitted, use standard output.
        type: string
      - id: reference
        doc: Reference file. Provide this to detect alleles through re-alignment.
            If no index (.fai) exists, it will be created.
        type:
          - File
          - 'null'
      - id: referenceFastaIndex
        doc: Index for the reference file.
        type:
          - File
          - 'null'
      - id: regions
        doc: Specify region(s) of interest to limit the tagging to reads/variants
            overlapping those regions. You can specify a space-separated list of regions
            in the form of chrom:start-end, chrom (consider entire chromosome), or
            chrom:start (consider region from this start to end of chromosome).
        type:
          - string
          - 'null'
      - id: sample
        doc: Name of a sample to phase. If not given, all samples in the input VCF
            are phased. Can be used multiple times.
        type:
          - string
          - 'null'
      - id: memory
        doc: The amount of memory this job will use.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 120
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: bam
        doc: BAM file containing tagged reads for haplotype.
        type: File
        outputBinding:
            glob: $(inputs.outputFile)
      - id: bamIndex
        doc: Index of the tagged BAM file.
        type: File
        outputBinding:
            glob: $(inputs.outputFile + ".bai")
