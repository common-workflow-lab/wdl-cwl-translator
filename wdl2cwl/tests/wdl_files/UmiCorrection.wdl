version 1.0

# Source: https://github.com/broadinstitute/warp/blob/11e9f8e6a4288dd3d8b23f60f7785a67ef1cbe31/tasks/skylab/UmiCorrection.wdl
#
# BSD 3-Clause License
# Copyright (c) 2021, Broad Institute
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

task CorrectUMItools {
    input {
        File bam_input
        File bam_index

        # runtime values
        String docker = "quay.io/humancellatlas/secondary-analysis-umitools:0.0.1"

        String output_bam_filename = "output.bam"
        String groupout_filename = "groupout.tsv"

        ## TODO: Optimize these values
        Int machine_mem_mb = 16000
        Int cpu = 1
        #Int disk = ceil(size(bam_input, "Gi") * 6) + 50
        Int preemptible = 3
    }

    meta {
        description: "Marks duplicates using umitools group specifically for single-cell experiments"
    }

    parameter_meta {
        bam_input: "Aligned and sorted bam"
        docker: "(optional) the docker image containing the runtime environment for this task"
        machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        #disk: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command {
        set -e

        mv ${bam_input} input.bam
        mv ${bam_index} input.bam.bai

        touch input.bam
        touch input.bam.bai

        umi_tools group \
            -I input.bam \
            -L outlog.txt \
            -E outerr.txt \
            -S duplicate_marked.bam \
            --output-bam \
            --extract-umi-method=tag \
            --umi-tag UR \
            --method directional \
            --per-gene \
            --per-cell \
            --cell-tag CB \
            --gene-tag GE \
            --no-sort-output \
            --group-out ${groupout_filename} \
            --umi-group-tag UB

       getUntaggedReads --in-bam-file input.bam --out-bam-file untagged.bam

       rm input.bam input.bam.bai
       samtools cat -o ${output_bam_filename} duplicate_marked.bam untagged.bam

    }

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        #disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File bam_output = "${output_bam_filename}"
        File group_output = "${groupout_filename}"
    }

}
