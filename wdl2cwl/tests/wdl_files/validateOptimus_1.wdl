version 1.0

# Source: https://github.com/broadinstitute/warp/blob/cec97750e3819fd88ba382534aaede8e05ec52df/tests/skylab/optimus/pr/ValidateOptimus.wdl
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

task ValidateBam {
    input {
        File bam
        String expected_checksum
    }

    Int required_disk = ceil(size(bam, "G") * 1.1)

    command <<<
        cacheInvalidationRandomString=4

        echo Starting checksum generation...

        # calculate hash for alignment positions only (a reduced bam hash)
        calculated_checksum=$( samtools view -F 256 "~{bam}" | cut -f 1-11 | md5sum | awk '{print $1}' )
        echo Reduced checksum generation complete

        if [ "$calculated_checksum" == "~{expected_checksum}" ]
        then
             echo Computed and expected bam hashes match \( "$calculated_checksum" \)
             printf PASS > result.txt
        else
             echo Computed \( "$calculated_checksum" \) and expected \( "~{expected_checksum}" \) bam file hashes do not match
             printf FAIL > result.txt
        fi
    >>>

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string("result.txt")
    }
}