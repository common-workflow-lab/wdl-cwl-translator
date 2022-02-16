cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: Format
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p \$(dirname $(inputs.outputPath))
                rtg RTG_MEM=$(inputs.rtgMem) format -f $(inputs.format) \
                -o $(inputs.outputPath) \
                $(inputs.inputFiles.map(function(el) {return el.path}).join(" "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/rtg-tools:3.10.1--0
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
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.inputFiles.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1 * 2)  * 60)
    inputs:
      - id: inputFiles
        doc: Input sequence files. May be specified 1 or more times.
        type:
            items: File
            type: array
      - id: format
        doc: Format of input. Allowed values are [fasta, fastq, fastq-interleaved,
            sam-se, sam-pe].
        default: fasta
        type: string
      - id: outputPath
        doc: Where the output should be placed.
        default: seq_data.sdf
        type: string
      - id: rtgMem
        doc: The amount of memory rtg will allocate to the JVM.
        default: 8G
        type: string
      - id: memory
        doc: The amount of memory this job will use.
        default: 9G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/rtg-tools:3.10.1--0
        type: string
    baseCommand:
      - bash
      - script.bash
    arguments:
      - valueFrom: ${if (inputs.inputFiles.length == 0) {throw "inputFiles must contain
            at least one item.";} else { return "";}}
    outputs:
      - id: sdf
        doc: RTGSequence Data File (SDF) format version of the input file(s).
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
  - cwlVersion: v1.2
    id: VcfEval
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputDir))"
                rtg RTG_MEM=$(inputs.rtgMem) vcfeval \
                --baseline $(inputs.baseline.path) \
                --calls $(inputs.calls.path) \
                $(inputs.evaluationRegions === null ? "" : "--evaluation-regions " + inputs.evaluationRegions.path) \
                $(inputs.bedRegions === null ? "" : "--bed-regions " + inputs.bedRegions.path) \
                --output $(inputs.outputDir) \
                --template $(inputs.template.path) \
                $(inputs.allRecords ? "--all-records" : "") \
                $(inputs.decompose ? "--decompose" : "") \
                $(inputs.refOverlap ? "--ref-overlap" : "") \
                $(inputs.sample === null ? "" : "--sample " + inputs.sample) \
                $(inputs.squashPloidy ? "--squash-ploidy" : "") \
                $("--output-mode " + inputs.outputMode) \
                --threads $(inputs.threads)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/rtg-tools:3.10.1--0
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
        timelimit: $(1 + Math.ceil((function(size_of=0){[inputs.baseline.path, inputs.calls.path].forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 5)  * 60)
    inputs:
      - id: baseline
        doc: VCF file containing baseline variants.
        type: File
      - id: baselineIndex
        doc: The baseline's VCF index.
        type: File
      - id: calls
        doc: VCF file containing called variants.
        type: File
      - id: callsIndex
        doc: The call's VCF index.
        type: File
      - id: squashPloidy
        doc: treat heterozygous genotypes as homozygous ALT in both baseline and calls,
            to allow matches that ignore zygosity differences.
        default: false
        type: boolean
      - id: outputMode
        doc: output reporting mode. Allowed values are [split, annotate, combine,
            ga4gh, roc-only] (Default is split).
        default: split
        type: string
      - id: outputDir
        doc: Directory for output.
        default: output/
        type: string
      - id: template
        doc: SDF of the reference genome the variants are called against.
        type: File
      - id: allRecords
        doc: use all records regardless of FILTER status (Default is to only process
            records where FILTER is "." or "PASS").
        default: false
        type: boolean
      - id: decompose
        doc: decompose complex variants into smaller constituents to allow partial
            credit.
        default: false
        type: boolean
      - id: refOverlap
        doc: allow alleles to overlap where bases of either allele are same-as-ref
            (Default is to only allow VCF anchor base overlap).
        default: false
        type: boolean
      - id: evaluationRegions
        doc: if set, evaluate within regions contained in the supplied BED file, allowing
            transborder matches. To be used for truth-set high-confidence regions
            or other regions of interest where region boundary effects should be minimized.
        type:
          - File
          - 'null'
      - id: bedRegions
        doc: if set, only read VCF records that overlap the ranges contained in the
            specified BED file.
        type:
          - File
          - 'null'
      - id: sample
        doc: the name of the sample to select. Use <baseline_sample>,<calls_sample>
            to select different sample names for baseline and calls. (Required when
            using multi-sample VCF files).
        type:
          - string
          - 'null'
      - id: rtgMem
        doc: The amount of memory rtg will allocate to the JVM.
        default: 8G
        type: string
      - id: threads
        doc: Number of threads. Default is 1.
        default: 1
        type: int
      - id: memory
        doc: The amount of memory this job will use.
        default: 9G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        type:
          - int
          - 'null'
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/rtg-tools:3.10.1--0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: falseNegativesVcf
        doc: Variants from thebaselineVCF which were not correctly called.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fn.vcf.gz")
      - id: falseNegativesVcfIndex
        doc: Index of the output VCF file `falseNegativesVcf`.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fn.vcf.gz.tbi")
      - id: falsePositivesVcf
        doc: Variants from thecallsVCF which do not agree with baseline variants.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fp.vcf.gz")
      - id: falsePositivesVcfIndex
        doc: Index of the output VCF file `falsePositivesVcf`.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fp.vcf.gz.tbi")
      - id: summary
        doc: Summary statistic file.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/summary.txt")
      - id: truePositivesBaselineVcf
        doc: Variants from thebaselineVCF which agree with variants in thecalls VCF.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp-baseline.vcf.gz")
      - id: truePositivesBaselineVcfIndex
        doc: Index of the output VCF file `truePositivesBaselineVcf`.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp-baseline.vcf.gz.tbi")
      - id: truePositivesVcf
        doc: Variants from thecallsVCF which agree with variants in the baseline VCF.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp.vcf.gz")
      - id: truePositivesVcfIndex
        doc: Index of the output VCF file `truePositivesVcf`.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp.vcf.gz.tbi")
      - id: nonSnpRoc
        doc: ROC data derived from those variants which were not represented asSNPs.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/non_snp_roc.tsv.gz")
      - id: phasing
        doc: Phasing file.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/phasing.txt")
      - id: weightedRoc
        doc: ROC data derived from all analyzed call variants, regardless of their
            representation.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/weighted_roc.tsv.gz")
      - id: allStats
        doc: All output files combined in a array.
        type:
            items: File
            type: array
        outputBinding:
            glob:
              - $(inputs.outputDir + "/fn.vcf.gz")
              - $(inputs.outputDir + "/fn.vcf.gz.tbi")
              - $(inputs.outputDir + "/fp.vcf.gz")
              - $(inputs.outputDir + "/fp.vcf.gz.tbi")
              - $(inputs.outputDir + "/tp-baseline.vcf.gz")
              - $(inputs.outputDir + "/tp-baseline.vcf.gz.tbi")
              - $(inputs.outputDir + "/tp.vcf.gz")
              - $(inputs.outputDir + "/tp.vcf.gz.tbi")
              - $(inputs.outputDir + "/summary.txt")
              - $(inputs.outputDir + "/non_snp_roc.tsv.gz")
              - $(inputs.outputDir + "/phasing.txt")
              - $(inputs.outputDir + "/weighted_roc.tsv.gz")
