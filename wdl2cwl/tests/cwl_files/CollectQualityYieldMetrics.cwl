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
              INPUT=$(inputs["input_bam"]===null?"":inputs["input_bam"]["class"]==="File"? "$(inputs.input_bam.path)": "$(inputs.input_bam)") \
              OQ=true \
              OUTPUT=$(inputs["metrics_filename"]===null?"":inputs["metrics_filename"]["class"]==="File"? "$(inputs.metrics_filename.path)": "$(inputs.metrics_filename)")
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 3584
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
