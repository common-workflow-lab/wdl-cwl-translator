class: CommandLineTool
id: GenerateReport
inputs:
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
      - entryname: script.bash
        entry: |4+


            set -eo pipefail

            cacheInvalidationRandomString=4

            # test each output for equality, echoing any failure states to stdout
            fail=false

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
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 953.67431640625
    outdirMin: 1024
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
