version 1.0

# Source: https://github.com/biowdl/tasks/blob/bc1bacf11498d2d30b85591cfccdcf71ef0966a5/rtg.wdl#L23
#
# Copyright (c) 2020 Leiden University Medical Center
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

task Format {
    input {
        Array[File]+ inputFiles
        String format = "fasta"
        String outputPath = "seq_data.sdf"

        String rtgMem = "8G"
        String memory = "9G"
        #Int timeMinutes = 1 + ceil(size(inputFiles) * 2)
        String dockerImage = "quay.io/biocontainers/rtg-tools:3.10.1--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        rtg RTG_MEM=~{rtgMem} format -f ~{format} \
        -o ~{outputPath} \
        ~{sep=' ' inputFiles}
    }

    output {
        File sdf = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFiles: {description: "Input sequence files. May be specified 1 or more times.", category: "required"}
        format: {description: "Format of input. Allowed values are [fasta, fastq, fastq-interleaved, sam-se, sam-pe].", category: "advanced"}
        outputPath: {description: "Where the output should be placed.", category: "advanced"}
        rtgMem: {description: "The amount of memory rtg will allocate to the JVM.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        sdf: {description: "RTGSequence Data File (SDF) format version of the input file(s)."}
    }
}
