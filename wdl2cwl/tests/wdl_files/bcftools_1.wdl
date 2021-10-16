version 1.0

# Copyright (c) 2018 Leiden University Medical Center
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

task Annotate {
    input {
        Array[String] columns = []
        Boolean force = false
        Boolean keepSites = false
        Boolean noVersion = false
        Array[String] samples = []
        Boolean singleOverlaps = false
        Array[String] removeAnns = []
        File inputFile
        String outputPath = "output.vcf.gz"

        File? annsFile
        String? collapse
        String? exclude
        File? headerLines
        String? newId
        String? include
        String? markSites
        String? regions
        File? regionsFile
        File? renameChrs
        File? samplesFile

        Int threads = 0
        String memory = "256M"
        # Int timeMinutes = 1 + ceil(size(inputFile, "G"))
        String dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    Boolean compressed = basename(outputPath) != basename(outputPath, ".gz")

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        bcftools annotate \
        -o ~{outputPath} \
        -O ~{true="z" false="v" compressed} \
        ~{"--annotations " + annsFile} \
        ~{"--collapse " + collapse} \
        ~{true="--columns" false="" length(columns) > 0} ~{sep="," columns} \
        ~{"--exclude " + exclude} \
        ~{true="--force" false="" force} \
        ~{"--header-lines " + headerLines} \
        ~{"--set-id " + newId} \
        ~{"--include " + include} \
        ~{true="--keep-sites" false="" keepSites} \
        ~{"--mark-sites " + markSites} \
        ~{true="--no-version" false="" noVersion} \
        ~{"--regions " + regions} \
        ~{"--regions-file " + regionsFile} \
        ~{"--rename-chrs " + renameChrs} \
        ~{true="--samples" false="" length(samples) > 0} ~{sep="," samples} \
        ~{"--samples-file " + samplesFile} \
        ~{true="--single-overlaps" false="" singleOverlaps} \
        ~{true="--remove" false="" length(removeAnns) > 0} ~{sep="," removeAnns} \
        ~{inputFile}

        ~{if compressed then 'bcftools index --tbi ~{outputPath}' else ''}
    }

    output {
        File outputVcf = outputPath
        File? outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        # time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        columns: {description: "Comma-separated list of columns or tags to carry over from the annotation file (see man page for details).", category: "advanced"}
        force: {description: "Continue even when parsing errors, such as undefined tags, are encountered.", category: "advanced"}
        keepSites: {description: "Keep sites which do not pass -i and -e expressions instead of discarding them.", category: "advanced"}
        noVersion: {description: "Do not append version and command line information to the output VCF header.", category: "advanced"}
        samples: {description: "List of samples for sample stats, \"-\" to include all samples.", category: "advanced"}
        singleOverlaps: {description: "keep memory requirements low with very large annotation files.", category: "advanced"}
        removeAnns: {description: "List of annotations to remove (see man page for details).", category: "advanced"}
        inputFile: {description: "A vcf or bcf file.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        annsFile: {description: "Bgzip-compressed and tabix-indexed file with annotations (see man page for details).", category: "advanced"}
        collapse: {description: "Treat as identical records with <snps|indels|both|all|some|none>, see man page for details.", category: "advanced"}
        exclude: {description: "Exclude sites for which the expression is true (see man page for details).", category: "advanced"}
        headerLines: {description: "Lines to append to the VCF header (see man page for details).", category: "advanced"}
        newId: {description: "Assign ID on the fly (e.g. --set-id +'%CHROM\_%POS').", category: "advanced"}
        include: {description: "Select sites for which the expression is true (see man page for details).", category: "advanced"}
        markSites: {description: "Annotate sites which are present ('+') or absent ('-') in the -a file with a new INFO/TAG flag.", category: "advanced"}
        regions: {description: "Restrict to comma-separated list of regions.", category: "advanced"}
        regionsFile: {description: "Restrict to regions listed in a file.", category: "advanced"}
        renameChrs: {description: "rename chromosomes according to the map in file (see man page for details).", category: "advanced"}
        samplesFile: {description: "File of samples to include.", category: "advanced"}
        threads: {description: "Number of extra decompression threads [0].", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVcf: {description: "Annotated VCF file."}
        outputVcfIndex: {description: "Index of the annotated VCF file."}
    }
}