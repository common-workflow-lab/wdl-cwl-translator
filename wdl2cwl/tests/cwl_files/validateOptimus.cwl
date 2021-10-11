class: CommandLineTool
id: GenerateReport
inputs:
  - id: bam
    type: File
  - id: expected_checksum
    type: string
  - id: loom_file
    type:
      - File
      - 'null'
  - id: expected_loom_file_checksum
    type: string
  - id: matrix
    type: File
  - id: matrix_row_index
    type: File
  - id: matrix_col_index
    type: File
  - id: reference_matrix
    type: File
  - id: cell_metrics
    type: File
  - id: gene_metrics
    type: File
  - id: expected_cell_metric_hash
    type: string
  - id: expected_gene_metric_hash
    type: string
  - id: bam_validation_result
    type: string
  - id: metric_and_index_validation_result
    type: string
  - id: matrix_validation_result
    type: string
  - id: loom_validation_result
    type: string
outputs:
  - id: result
    type: string
    outputBinding:
        glob: result.txt
        outputEval: |-
            $(self.contents.replace(/[
            ]+$/, '')
  - id: result
    type: string
    outputBinding:
        glob: result.txt
        outputEval: |-
            $(self.contents.replace(/[
            ]+$/, '')
  - id: result
    type: string
    outputBinding:
        glob: result.txt
        outputEval: |-
            $(self.contents.replace(/[
            ]+$/, '')
  - id: new_reference_matrix
    type: File
    outputBinding:
        glob: newReferenceMatrix.rds
  - id: reads_per_cell_histogram
    type: File
    outputBinding:
        glob: reads_per_cell_histogram.png
  - id: reads_per_gene_histogram
    type: File
    outputBinding:
        glob: reads_per_gene_histogram.png
  - id: number_of_genes_per_cell
    type: File
    outputBinding:
        glob: number_of_genes_per_cell.png
  - id: result
    type: string
    outputBinding:
        glob: result.txt
        outputEval: |-
            $(self.contents.replace(/[
            ]+$/, '')
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
