cwlVersion: v1.2
$graph:
  - class: CommandLineTool
    id: Format
    inputs:
      - id: inputFiles
        type:
            items: File
            type: array
      - id: format
        default: fasta
        type: string
      - id: outputPath
        default: seq_data.sdf
        type: string
      - id: rtgMem
        default: 8G
        type: string
      - id: memory
        default: 9G
        type: string
      - id: timeMinutes
        type:
          - int
          - 'null'
      - id: dockerImage
        default: quay.io/biocontainers/rtg-tools:3.10.1--0
        type: string
    outputs:
      - id: sdf
        type: File
        outputBinding:
            glob: $(inputs.outputPath)
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){inputs.inputFiles.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1 * 2)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: VcfEval
    inputs:
      - id: baseline
        type: File
      - id: baselineIndex
        type: File
      - id: calls
        type: File
      - id: callsIndex
        type: File
      - id: squashPloidy
        default: false
        type: boolean
      - id: outputMode
        default: split
        type: string
      - id: outputDir
        default: output/
        type: string
      - id: template
        type: File
      - id: allRecords
        default: false
        type: boolean
      - id: decompose
        default: false
        type: boolean
      - id: refOverlap
        default: false
        type: boolean
      - id: evaluationRegions
        type:
          - File
          - 'null'
      - id: bedRegions
        type:
          - File
          - 'null'
      - id: sample
        type:
          - string
          - 'null'
      - id: rtgMem
        default: 8G
        type: string
      - id: threads
        default: 1
        type: int
      - id: memory
        default: 9G
        type: string
      - id: timeMinutes
        type:
          - int
          - 'null'
      - id: dockerImage
        default: quay.io/biocontainers/rtg-tools:3.10.1--0
        type: string
    outputs:
      - id: falseNegativesVcf
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fn.vcf.gz")
      - id: falseNegativesVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fn.vcf.gz.tbi")
      - id: falsePositivesVcf
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fp.vcf.gz")
      - id: falsePositivesVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/fp.vcf.gz.tbi")
      - id: summary
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/summary.txt")
      - id: truePositivesBaselineVcf
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp-baseline.vcf.gz")
      - id: truePositivesBaselineVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp-baseline.vcf.gz.tbi")
      - id: truePositivesVcf
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp.vcf.gz")
      - id: truePositivesVcfIndex
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/tp.vcf.gz.tbi")
      - id: nonSnpRoc
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/non_snp_roc.tsv.gz")
      - id: phasing
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/phasing.txt")
      - id: weightedRoc
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "/weighted_roc.tsv.gz")
      - id: allStats
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
                --output-mode  $(inputs.outputMode) \
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(1 + Math.ceil((function(size_of=0){[inputs.baseline.path, inputs.calls.path].forEach(function(element){
            if (element) {size_of += element.size}})}) / 1000^3 * 5)  * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
