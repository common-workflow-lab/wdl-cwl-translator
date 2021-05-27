#!/usr/bin/env cwl-runner

class: CommandLineTool
id: CollectQualityYieldMetrics
cwlVersion: v1.0

requirements:
- class: DockerRequirement
  dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
- class: InlineJavascriptRequirement

hints:
- class: ResourceRequirement
  ramMin: 3584
  
inputs:
  input_bam:
    type: File
    inputBinding:
      position: 1
  
  metrics_filename:
    type: string
    inputBinding:
      position: 2

outputs:
  quality_yield_metrics:
    type: File
    outputBinding:
      glob: $(inputs.metrics_filename)

baseCommand: [java -Xms2000m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=$(input_bam) \
      OQ=true \
      OUTPUT=$(metrics_filename)]