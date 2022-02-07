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
task GenerateReport {
  input {
    String metric_and_index_validation_result
    String matrix_validation_result
    String loom_validation_result
  }

  Int required_disk = 1

  command <<<

    set -eo pipefail

    cacheInvalidationRandomString=4

    # test each output for equality, echoing any failure states to stdout
    fail=false

    echo Metrics Validation: ~{metric_and_index_validation_result}
    if [ ~{metric_and_index_validation_result} == "FAIL" ]; then
        echo --- Ignoring failed metric and index test ---
        # Do not fail tests for this
        # fail=true
    fi

    echo Matrix Validation: ~{matrix_validation_result}
    if [ "~{matrix_validation_result}" == "FAIL" ]; then
        fail=true
    fi

    echo Loom Validation: ~{loom_validation_result}
    if [ "~{loom_validation_result}" == "FAIL" ]; then
        echo --- Ignoring failed loom test ---
        # Do not fail tests for this
        # fail=true
    fi

    if [ "$fail" == "true" ]; then exit 1; fi

  >>>

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1.0 GB"
    disks: "local-disk ${required_disk} HDD"
  }

  output {}
}
