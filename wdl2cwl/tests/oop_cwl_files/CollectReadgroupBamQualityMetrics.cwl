class: CommandLineTool
id: CollectReadgroupBamQualityMetrics
inputs:
  - id: input_bam
    type: File
  - id: input_bam_index
    type: File
  - id: output_bam_prefix
    type: string
  - id: ref_dict
    type: File
  - id: ref_fasta
    type: File
  - id: ref_fasta_index
    type: File
  - id: collect_gc_bias_metrics
    default: true
    type: boolean
  - id: preemptible_tries
    type: int
outputs:
  - id: alignment_summary_metrics
    type: File
    outputBinding:
        glob: $(inputs.output_bam_prefix + '.alignment_summary_metrics')
  - id: gc_bias_detail_metrics
    type: File
    outputBinding:
        glob: $(inputs.output_bam_prefix + '.gc_bias.detail_metrics')
  - id: gc_bias_pdf
    type: File
    outputBinding:
        glob: $(inputs.output_bam_prefix + '.gc_bias.pdf')
  - id: gc_bias_summary_metrics
    type: File
    outputBinding:
        glob: $(inputs.output_bam_prefix + '.gc_bias.summary_metrics')
requirements:
  - class: DockerRequirement
    dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            # These are optionally generated, but need to exist for Cromwell's sake
            touch $(inputs.output_bam_prefix).gc_bias.detail_metrics \
            $(inputs.output_bam_prefix).gc_bias.pdf \
            $(inputs.output_bam_prefix).gc_bias.summary_metrics

            java -Xms5000m -jar /usr/picard/picard.jar \
              CollectMultipleMetrics \
              INPUT=$(inputs.input_bam.path) \
            REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
            OUTPUT=$(inputs.output_bam_prefix) \
            ASSUME_SORTED=true \
            PROGRAM=null \
            PROGRAM=CollectAlignmentSummaryMetrics \
            $(inputs.collect_gc_bias_metrics ? "PROGRAM="CollectGcBiasMetrics"" : "") \
              METRIC_ACCUMULATION_LEVEL=null \
              METRIC_ACCUMULATION_LEVEL=READ_GROUP
            sed -i -e 1,5d "$(inputs.output_bam_prefix).alignment_summary_metrics"   # for reproducibility
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    ramMin: 7168.0
    outdirMin: 30
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
