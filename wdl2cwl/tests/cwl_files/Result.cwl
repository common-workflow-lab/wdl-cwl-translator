#!/usr/bin/env cwl-runner

class: CommandLineTool
id: CollectQualityYieldMetrics
cwlVersion: v1.0

baseCommand: ["sh", "example.sh"]

requirements:
- class: InitialWorkDirRequirement
  listing:
  - entryname: example.sh
    entry: 'java -Xms2000m -jar /usr/picard/picard.jar \    CollectQualityYieldMetrics \    INPUT=$(inputs.input_bam.path) \    OQ=true \    OUTPUT=$(inputs.metrics_filename)'
- class: DockerRequirement
  dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
- class: InlineJavascriptRequirement

hints:
- class: ResourceRequirement
  ramMin: 3584
  
inputs:
  input_bam:
    type: File
  
  metrics_filename:
    type: string

outputs:
  quality_yield_metrics:
    type: File
    outputBinding:
      glob: $(inputs.metrics_filename)


