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

task ValidateMetrics {
    input {
        File cell_metrics
        File gene_metrics

        String expected_cell_metric_hash
        String expected_gene_metric_hash
    }

    Int required_disk = ceil( (size(cell_metrics, "G") + size(gene_metrics, "G") )* 1.1)

    command <<<
        set -eo pipefail

        cacheInvalidationRandomString=4

        # check matrix row and column indexes files hash
        gene_metric_hash=$(zcat "~{gene_metrics}" | md5sum | awk '{print $1}')
        cell_metric_hash=$(zcat "~{cell_metrics}" | md5sum | awk '{print $1}')

        fail=false

        if [ "$gene_metric_hash" == "~{expected_gene_metric_hash}" ]; then
            echo Computed and expected gene metrics match \( "$gene_metric_hash" \)
        else
            echo Computed \( "$gene_metric_hash" \) and expected \( "~{expected_gene_metric_hash}" \) gene checksums do not match
            fail=true
        fi

        if [ "$cell_metric_hash" == "~{expected_cell_metric_hash}" ]; then
            echo Computed and expected cell metrics match \( "$cell_metric_hash" \)
        else
            echo Computed \( "$cell_metric_hash" \) and expected \( "~{expected_cell_metric_hash}" \) cell metrics hashes do not match
            fail=true
        fi

        if [ $fail == "true" ]; then
            printf FAIL > result.txt
        else
            printf PASS > result.txt
        fi
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "1.00 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string("result.txt")
    }

}
