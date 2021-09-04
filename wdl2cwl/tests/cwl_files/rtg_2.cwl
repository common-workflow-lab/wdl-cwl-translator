class: CommandLineTool
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
  - id: template
    type: File
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
  - id: squashPloidy
    default: false
    type: boolean
  - id: outputMode
    default: split
    type: string
  - id: outputDir
    default: output/
    type: string
  - id: allRecords
    default: false
    type: boolean
  - id: decompose
    default: false
    type: boolean
  - id: refOverlap
    default: false
    type: boolean
  - id: rtgMem
    default: 8G
    type: string
  - id: threads
    default: 1
    type: int
  - id: memory
    default: 9G
    type: string
  - id: dockerImage
    default: quay.io/biocontainers/rtg-tools:3.10.1--0
    type: string
outputs:
  - id: falseNegativesVcf
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/fn.vcf.gz
  - id: falseNegativesVcfIndex
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/fn.vcf.gz.tbi
  - id: falsePositivesVcf
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/fp.vcf.gz
  - id: falsePositivesVcfIndex
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/fp.vcf.gz.tbi
  - id: summary
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/summary.txt
  - id: truePositivesBaselineVcf
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/tp-baseline.vcf.gz
  - id: truePositivesBaselineVcfIndex
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/tp-baseline.vcf.gz.tbi
  - id: truePositivesVcf
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/tp.vcf.gz
  - id: truePositivesVcfIndex
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/tp.vcf.gz.tbi
  - id: nonSnpRoc
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/non_snp_roc.tsv.gz
  - id: phasing
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/phasing.txt
  - id: weightedRoc
    type: File
    outputBinding:
        glob: $(inputs.outputDir)/weighted_roc.tsv.gz
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/rtg-tools:3.10.1--0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p "\$(dirname $(inputs.outputDir))"
            rtg RTG_MEM=$(inputs.rtgMem) vcfeval \
            --baseline $(inputs.baseline.path) \
            --calls $(inputs.calls.path) \
            $(inputs["evaluationRegions"].path === null ? "" : "--evaluation-regions inputs["evaluationRegions"].path") \
            $(inputs["bedRegions"].path === null ? "" : "--bed-regions inputs["bedRegions"].path") \
            --output $(inputs.outputDir) \
            --template $(inputs.template.path) \
            $(inputs["allRecords"] ? "--all-records" : "") \
            $(inputs["decompose"] ? "--decompose" : "") \
            $(inputs["refOverlap"] ? "--ref-overlap" : "") \
            $(inputs["sample"] === null ? "" : "--sample inputs["sample"]") \
            $(inputs["squashPloidy"] ? "--squash-ploidy" : "") \
            --output-mode $(inputs.outputMode) \
            --threads $(inputs.threads)
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: |-
        ${
        var unit = inputs["memory"].match(/[a-zA-Z]+/g).join("");
        var value = parseInt(inputs["memory"].match(/[0-9]+/g));
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
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
