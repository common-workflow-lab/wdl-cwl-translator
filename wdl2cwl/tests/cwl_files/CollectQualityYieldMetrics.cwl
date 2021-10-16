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
        entry: |4

            java -Xms2000m -jar /usr/picard/picard.jar \
              CollectQualityYieldMetrics \
              INPUT=$(inputs.input_bam.path) \
              OQ=true \
              OUTPUT=$(inputs.metrics_filename)
            sed -i -e 1,5d $(inputs.metrics_filename)  # for reproducibility
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 3584
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
