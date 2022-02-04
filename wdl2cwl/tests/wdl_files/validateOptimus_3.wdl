version 1.0

# Source: https://github.com/broadinstitute/warp/blob/24274db3c6de25b1cdaecf4bc2f8d16be554d3e8/tests/skylab/optimus/pr/ValidateOptimus.wdl
#
# Copyright Broad Institute, 2020
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task ValidateMatrix {
    input {
        File matrix
        File matrix_row_index
        File matrix_col_index
        File reference_matrix
    }

    Int required_disk = ceil( size(matrix, "G") * 1.1 )

    command <<<
        cacheInvalidationRandomString=4
       
       ## Convert matrix to format that can be read by R
       npz2rds.sh -c ~{matrix_col_index} -r ~{matrix_row_index} \
           -d ~{matrix} -o matrix.rds

       cp ~{reference_matrix} referenceMatrix.rds

       ## Run tests
       Rscript /root/tools/checkMatrix.R
       checkMatrixResult=$?

       if [ $checkMatrixResult == 0 ]; then
           printf PASS > result.txt
       else
           printf FAIL > result.txt
       fi

    >>>

    runtime {
        docker: "quay.io/humancellatlas/optimus-matrix-test:0.0.7"
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string('result.txt')
        File new_reference_matrix = "newReferenceMatrix.rds"
        File reads_per_cell_histogram = "reads_per_cell_histogram.png"
        File reads_per_gene_histogram = "reads_per_gene_histogram.png"
        File number_of_genes_per_cell = "number_of_genes_per_cell.png"
    }

}
