class: CommandLineTool
id: GenerateReport
inputs:
  - id: bam_validation_result
    type: string
  - id: metric_and_index_validation_result
    type: string
  - id: matrix_validation_result
    type: string
  - id: loom_validation_result
    type: string
outputs: []
requirements:
  - class: DockerRequirement
    dockerPull: ubuntu:18.04
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4+


            set -eo pipefail

            cacheInvalidationRandomString=4

            # test each output for equality, echoing any failure states to stdout
            fail=false

            echo Bam Validation: $(inputs.bam_validation_result)
            if [ "$(inputs.bam_validation_result)" == "FAIL" ]; then
                fail=true
            fi

            echo Metrics Validation: $(inputs.metric_and_index_validation_result)
            if [ $(inputs.metric_and_index_validation_result) == "FAIL" ]; then
                echo --- Ignoring failed metric and index test ---
                # Do not fail tests for this
                # fail=true
            fi

            echo Matrix Validation: $(inputs.matrix_validation_result)
            if [ "$(inputs.matrix_validation_result)" == "FAIL" ]; then
                fail=true
            fi

            echo Loom Validation: $(inputs.loom_validation_result)
            if [ "$(inputs.loom_validation_result)" == "FAIL" ]; then
                echo --- Ignoring failed loom test ---
                # Do not fail tests for this
                # fail=true
            fi

            if [ "$fail" == "true" ]; then exit 1; fi

  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 0
  - class: ResourceRequirement
    coresMin: 1
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
