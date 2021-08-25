class: CommandLineTool
id: RunDeepVariant
inputs:
  - id: referenceFasta
    type: File
  - id: referenceFastaIndex
    type: File
  - id: inputBam
    type: File
  - id: inputBamIndex
    type: File
  - id: modelType
    type: string
  - id: outputVcf
    type: string
  - id: postprocessVariantsExtraArgs
    default: ''
    type:
      - string
      - 'null'
  - id: customizedModel
    default: ''
    type:
      - File
      - 'null'
  - id: numShards
    default: ''
    type:
      - int
      - 'null'
  - id: outputGVcf
    default: ''
    type:
      - string
      - 'null'
  - id: outputGVcfIndex
    default: ''
    type:
      - string
      - 'null'
  - id: regions
    default: ''
    type:
      - File
      - 'null'
  - id: sampleName
    default: ''
    type:
      - string
      - 'null'
  - id: VCFStatsReport
    default: true
    type:
      - boolean
      - 'null'
  - id: memory
    default: 3G
    type: string
  - id: timeMinutes
    default: 5000
    type: int
  - id: dockerImage
    default: google/deepvariant:1.0.0
    type: string
outputs:
  - id: outputVCF
    type: File
    outputBinding:
        glob: $(inputs.outputVcf)
  - id: outputVCFIndex
    type: File
    outputBinding:
        glob: $(inputs.outputVcf).tbi
  - id: outputVCFStatsReport
    type:
      - items: File
        type: array
    outputBinding:
        glob: '*.visual_report.html'
  - id: outputGVCF
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputGVcf)
  - id: outputGVCFIndex
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputGVcfIndex)
requirements:
  - class: DockerRequirement
    dockerPull: google/deepvariant:1.0.0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            /opt/deepvariant/bin/run_deepvariant \
            --ref $(inputs.referenceFasta.path) \
            --reads $(inputs.inputBam.path) \
            --model_type $(inputs.modelType) \
            --output_vcf $(inputs.outputVcf) \
            --output_gvcf $(inputs.outputGVcf) \
            --customized_model $(inputs.customizedModel) \
            --num_shards $(inputs.numShards) \
            --regions $(inputs.regions) \
            --sample_name $(inputs.sampleName) \
            --postprocess_variants_extra_args $(inputs.postprocessVariantsExtraArgs) \
            $(inputs["VCFStatsReport"] ? "--vcf_stats_report" : "--novcf_stats_report")
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
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes* 60)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
