cwlVersion: v1.2
id: RunDeepVariant
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            set -e
            mkdir reference_dir
            ln -s $(inputs.referenceFasta.path) reference_dir/$(inputs.referenceFasta.basename)
            ln -s $(inputs.referenceFastaIndex.path) reference_dir/$(inputs.referenceFastaIndex.basename)
            mkdir bam_dir
            ln -s $(inputs.inputBam.path) bam_dir/$(inputs.inputBam.basename)
            ln -s $(inputs.inputBamIndex.path) bam_dir/$(inputs.inputBamIndex.basename)
            /opt/deepvariant/bin/run_deepvariant \
            --ref reference_dir/$(inputs.referenceFasta.basename) \
            --reads bam_dir/$(inputs.inputBam.basename) \
            --model_type $(inputs.modelType) \
            --output_vcf $(inputs.outputVcf) \
            $(inputs.outputGVcf === null ? "" : "--output_gvcf " + inputs.outputGVcf) \
            $(inputs.customizedModel === null ? "" : "--customized_model " + inputs.customizedModel.path) \
            $(inputs.numShards === null ? "" : "--num_shards " + inputs.numShards) \
            $("--regions " + inputs.regions) \
            $(inputs.sampleName === null ? "" : "--sample_name " + inputs.sampleName) \
            $(inputs.postprocessVariantsExtraArgs === null ? "" : "--postprocess_variants_extra_args " + inputs.postprocessVariantsExtraArgs) \
            $(inputs.VCFStatsReport === null ? "" : inputs.VCFStatsReport ? "--vcf_stats_report" : "--novcf_stats_report")
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: google/deepvariant:1.0.0
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
  - id: referenceFasta
    doc: Genome reference to use.
    type: File
  - id: referenceFastaIndex
    doc: Index for the genome reference file.
    type: File
  - id: inputBam
    doc: Aligned, sorted, indexed BAM file containing the reads we want to call.
    type: File
  - id: inputBamIndex
    doc: Index for the input bam file.
    type: File
  - id: modelType
    doc: <WGS|WES|PACBIO>. Type of model to use for variant calling. Each model_type
        has an associated default model, which can be overridden by the --customized_model
        flag.
    type: string
  - id: outputVcf
    doc: Path where we should write VCF file.
    type: string
  - id: postprocessVariantsExtraArgs
    doc: A comma-separated list of flag_name=flag_value. 'flag_name' has to be valid
        flags for calpostprocess_variants.py.
    type:
      - string
      - 'null'
  - id: customizedModel
    doc: A path to a model checkpoint to load for the `call_variants` step. If not
        set, the default for each --model_type will be used.
    type:
      - File
      - 'null'
  - id: numShards
    doc: Number of shards for make_examples step.
    type:
      - int
      - 'null'
  - id: outputGVcf
    doc: Path where we should write gVCF file.
    type:
      - string
      - 'null'
  - id: outputGVcfIndex
    doc: Path to where the gVCF index file will be written. This is needed as a workaround,
        set it to `outputGVcf + '.tbi.'`
    type:
      - string
      - 'null'
  - id: regions
    doc: List of regions we want to process, in BED/BEDPE format.
    type: string
  - id: sampleName
    doc: Sample name to use instead of the sample name from the input reads BAM (SM
        tag in the header).
    type:
      - string
      - 'null'
  - id: VCFStatsReport
    doc: Output a visual report (HTML) of statistics about the output VCF.
    default: true
    type:
      - boolean
      - 'null'
  - id: memory
    doc: The amount of memory this job will use.
    default: 3G
    type: string
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    default: 5000
    type: int
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: google/deepvariant:1.0.0
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: outputVCF
    doc: Output VCF file.
    type: File
    outputBinding:
        glob: $(inputs.outputVcf)
  - id: outputVCFIndex
    doc: Index of output VCF file.
    type: File
    outputBinding:
        glob: $(inputs.outputVcf + ".tbi")
  - id: outputVCFStatsReport
    doc: Statistics file.
    type:
        items: File
        type: array
    outputBinding:
        glob: $("*.visual_report.html")
  - id: outputGVCF
    doc: GVCF version of VCF file(s).
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputGVcf)
  - id: outputGVCFIndex
    doc: Index of GVCF file(s).
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputGVcfIndex)
