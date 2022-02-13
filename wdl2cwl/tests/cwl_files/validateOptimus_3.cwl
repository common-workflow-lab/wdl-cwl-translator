cwlVersion: v1.2
id: ValidateMatrix
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
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
hints:
  - class: DockerRequirement
    dockerPull: quay.io/humancellatlas/optimus-matrix-test:0.0.7
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 15258.7890625
    outdirMin: $((Math.ceil((function(size_of=0){inputs.matrix.path.forEach(function(element){
        if (element) {size_of += element.size}})}) / 1000^3 * 1.1) ) * 1024)
inputs:
  - id: matrix
    type: File
  - id: matrix_row_index
    type: File
  - id: matrix_col_index
    type: File
  - id: reference_matrix
    type: File
baseCommand:
  - bash
  - script.bash
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
