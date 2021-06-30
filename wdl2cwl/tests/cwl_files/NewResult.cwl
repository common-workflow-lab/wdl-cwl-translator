class: CommandLineTool
id: CollectQualityYieldMetrics
inputs:
- id: input_bam
  type: File
- id: metrics_filename
  type: string
- id: preemptible_tries
  type: int
outputs:
- id: quality_yield_metrics
  type: File
  outputBinding:
    glob: $(inputs.metrics_filename)
requirements:
- class: DockerRequirement
  dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
- class: InitialWorkDirRequirement
  listing:
  - entryname: example.sh
    entry: "\njava -Xms2000m -jar /usr/picard/picard.jar \n  CollectQualityYieldMetrics\
      \ \n  INPUT=$(inputs.input_bam.path) \n  OQ=true \n  OUTPUT=$(inputs.metrics_filename)\n"
- class: InlineJavascriptRequirement
hints:
- class: ResourceRequirement
  ramMin: 3584
cwlVersion: v1.2
baseCommand:
- sh
- example.sh
