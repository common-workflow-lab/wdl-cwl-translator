class: CommandLineTool
id: ValidateMatrix
inputs:
  - id: matrix
    type: File
  - id: matrix_row_index
    type: File
  - id: matrix_col_index
    type: File
  - id: reference_matrix
    type: File
outputs:
  - id: result
    type: string
    outputBinding:
        loadContents: true
        glob: result.txt
        outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
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
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/humancellatlas/optimus-matrix-test:0.0.7
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4+

             cacheInvalidationRandomString=4

            ## Convert matrix to format that can be read by R
            npz2rds.sh -c $(inputs.matrix_col_index.path) -r $(inputs.matrix_row_index.path) \
            -d $(inputs.matrix.path) -o matrix.rds

            cp $(inputs.reference_matrix.path) referenceMatrix.rds

            ## Run tests
            Rscript /root/tools/checkMatrix.R
            checkMatrixResult=$?

            if [ $checkMatrixResult == 0 ]; then
                printf PASS > result.txt
            else
                printf FAIL > result.txt
            fi

  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 15258.7890625
    outdirMin: 1024
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
