class: CommandLineTool
id: CollectQualityYieldMetrics
inputs:
- id: input_bam
  type: File
- id: metrics_filename
  type: string
outputs:
- id: quality_yield_metrics
  type: File
  outputBinding:
    glob: $(inputs.metrics_filename)
requirements:
- class: DockerRequirement
  dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
hints:
- class: ResourceRequirement
  ramMin: 3584
cwlVersion: v1.0
baseCommand:
- java
- -Xms2000m
- -jar
- /usr/picard/picard.jar
- CollectQualityYieldMetrics
arguments:
- valueFrom: INPUT=$(inputs.input_bam.path)
- valueFrom: OQ=true
- valueFrom: OUTPUT=$(inputs.metrics_filename)

