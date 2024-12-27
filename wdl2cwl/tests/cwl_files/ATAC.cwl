cwlVersion: v1.2
id: ATAC
class: Workflow
doc: Processing for single-cell ATAC-seq data from the level of raw fastq reads to
    the generation of a snap file with snaptools. ATAC-seq (Assay for Transposase-Accessible
    Chromatin using sequencing) is a technique used in molecular biology to assess
    genome-wide chromatin accessibility. This pipeline accepts fastq files where the
    cell barcode has been added to the fastq read names as the first field.
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
inputs:
  - id: fastq_gzipped_input_read1
    doc: read 1 fastq file as input for the pipeline, the cellular barcodes must be
        the first part of the read name seperated by colon
    type: File
  - id: fastq_gzipped_input_read2
    doc: read 2 fastq file as input for the pipeline, the cellular barcodes must be
        the first part of the read name separated by colon
    type: File
  - id: min_length
    doc: minimum length for trimming. Reads that are too short even before adapter
        removal are also discarded
    type: int
  - id: quality_cutoff
    doc: cutadapt option to trim low-quality ends from reads before adapter removal
    type: int
  - id: adapter_seq_read1
    doc: cutadapt option for the sequence adapter for read 1 fastq
    type: string
  - id: adapter_seq_read2
    doc: cutadapt option for the sequence adapter for read 2 fastq
    type: string
  - id: tar_bwa_reference
    doc: the pre built tar file containing the reference fasta and corresponding reference
        files for the BWA aligner
    type: File
  - id: read_group_id
    doc: the read group id to be added upon alignment
    default: RG1
    type: string
  - id: read_group_sample_name
    doc: the read group sample to be added upon alignment
    default: RGSN1
    type: string
  - id: bwa_cpu
    doc: the number of cpu cores to use during alignment
    default: 16
    type: int
  - id: genome_name
    doc: the name of the genome being analyzed, input to snap tools, curently mm10
        and hg19 supported
    type: string
  - id: genome_size_file
    doc: name of the file with chromosome sizes for the genome in the input tar file
    type: File
  - id: min_map_quality
    doc: the minimum mapping quality to be filtered by samtools view and snap-pre
        (snaptools task)
    type: int
  - id: max_fragment_length
    doc: the maximum fragment length for filtering out reads by gatk and snap-pre
        (snaptools task)
    type: int
  - id: output_base_name
    doc: base name to be used for the pipelines output and intermediate files
    type: string
  - id: bin_size_list
    doc: space separated list of bins to generate
    default: '10000'
    type: string
  - id: TrimAdapters.docker_image
    default: quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1
    type: string
  - id: BWAPairedEndAlignment.docker_image
    default: quay.io/humancellatlas/snaptools:0.0.1
    type: string
  - id: SamToBam.docker_image
    default: quay.io/biocontainers/samtools:1.9--h10a08f8_12
    type: string
  - id: SortCoordinateOrder.sort_order
    default: coordinate
    type: string
  - id: SortCoordinateOrder.docker_image
    default: quay.io/biocontainers/picard:2.18.3--py27_0
    type: string
  - id: FilterMinMapQuality.docker_image
    default: quay.io/biocontainers/samtools:1.9--h10a08f8_12
    type: string
  - id: FilterMaxFragmentLength.docker_image
    default: broadinstitute/gatk:4.1.2.0
    type: string
  - id: FilterMitochondrialReads.docker_image
    default: quay.io/biocontainers/samtools:1.9--h10a08f8_12
    type: string
  - id: MakeCompliantChrMBAM.docker_image
    default: quay.io/humancellatlas/snaptools:0.0.1
    type: string
  - id: SortQueryName.docker_image
    default: quay.io/biocontainers/picard:2.18.3--py27_0
    type: string
  - id: MakeCompliantFilteredAndSortedBAM.docker_image
    default: quay.io/humancellatlas/snaptools:0.0.1
    type: string
  - id: SnapPre.docker_image
    default: quay.io/humancellatlas/snaptools:0.0.1
    type: string
  - id: SnapCellByBin.snap_output_name
    default: output.snap
    type: string
  - id: SnapCellByBin.docker_image
    default: quay.io/humancellatlas/snaptools:0.0.1
    type: string
  - id: BreakoutSnap.docker_image
    default: quay.io/humancellatlas/snap-breakout:0.0.1
    type: string
steps:
  - id: TrimAdapters
    in:
      - id: fastq_input_read1
        source: fastq_gzipped_input_read1
      - id: fastq_input_read2
        source: fastq_gzipped_input_read2
      - id: min_length
        source: min_length
      - id: quality_cutoff
        source: quality_cutoff
      - id: adapter_seq_read1
        source: adapter_seq_read1
      - id: adapter_seq_read2
        source: adapter_seq_read2
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: TrimAdapters.docker_image
    out:
      - id: fastq_trimmed_adapter_output_read1
      - id: fastq_trimmed_adapter_output_read2
      - id: monitoring_log
    run:
        id: TrimAdapters
        class: CommandLineTool
        inputs:
          - id: fastq_input_read1
            doc: read 1 fastq file as input for the pipeline
            type: File
          - id: fastq_input_read2
            doc: read 2 fastq file as input for the pipeline
            type: File
          - id: min_length
            doc: the minimum legnth for trimming. Reads that are too short even before
                adapter removal are also discarded
            type: int
          - id: quality_cutoff
            doc: cutadapt option to trim low-quality ends from reads before adapter
                removal
            type: int
          - id: adapter_seq_read1
            doc: cutadapt option for the sequence adapter for read 1 fastq
            type: string
          - id: adapter_seq_read2
            doc: cutadapt option for the sequence adapter for read 2 fastq
            type: string
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using cutadapt to be used (default: quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1)'
            default: quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1
            type: string
        outputs:
          - id: fastq_trimmed_adapter_output_read1
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".R1.trimmed_adapters.fastq.gz")
          - id: fastq_trimmed_adapter_output_read2
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".R2.trimmed_adapters.fastq.gz")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    # fastq's, "-f", -A for paired adapters read 2"
                    cutadapt \
                      -f fastq \
                      --minimum-length $(inputs.min_length) \
                      --quality-cutoff $(inputs.quality_cutoff) \
                      --adapter $(inputs.adapter_seq_read1) \
                      -A $(inputs.adapter_seq_read2) \
                      --output $(inputs.output_base_name + ".R1.trimmed_adapters.fastq.gz") \
                      --paired-output $(inputs.output_base_name + ".R2.trimmed_adapters.fastq.gz") \
                      $(inputs.fastq_input_read1.path) $(inputs.fastq_input_read2.path)
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3840.0
            outdirMin: '$((Math.ceil(2 * (function(size_of=0){inputs.fastq_input_read1.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.fastq_input_read2.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.fastq_input_read1.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.fastq_input_read2.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BWAPairedEndAlignment
    in:
      - id: fastq_input_read1
        source: TrimAdapters/fastq_trimmed_adapter_output_read1
      - id: fastq_input_read2
        source: TrimAdapters/fastq_trimmed_adapter_output_read2
      - id: tar_bwa_reference
        source: tar_bwa_reference
      - id: read_group_id
        source: read_group_id
      - id: read_group_sample_name
        source: read_group_sample_name
      - id: cpu
        source: bwa_cpu
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: BWAPairedEndAlignment.docker_image
    out:
      - id: sam_aligned_output
      - id: monitoring_log
    run:
        id: BWAPairedEndAlignment
        class: CommandLineTool
        inputs:
          - id: fastq_input_read1
            doc: the trimmed read 1 fastq file as input for the aligner
            type: File
          - id: fastq_input_read2
            doc: the trimmed read 1 fastq file as input for the aligner
            type: File
          - id: tar_bwa_reference
            doc: the pre built tar file containing the reference fasta and cooresponding
                reference files for the BWA aligner
            type: File
          - id: read_group_id
            doc: the read group id to be added upon alignment
            type: string
          - id: read_group_sample_name
            doc: the read group sample to be added upon alignment
            type: string
          - id: cpu
            doc: the number of cpu cores to use during alignment
            type: int
          - id: output_base_name
            doc: basename to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using BWA to be used (default: quay.io/humancellatlas/snaptools:0.0.1)'
            default: quay.io/humancellatlas/snaptools:0.0.1
            type: string
        outputs:
          - id: sam_aligned_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".aligned.sam")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    # prepare reference
                    declare -r REF_DIR=\$(mktemp -d genome_referenceXXXXXX)
                    tar -xf "$(inputs.tar_bwa_reference.path)" -C $REF_DIR --strip-components 1
                    rm "$(inputs.tar_bwa_reference.path)" || /bin/true

                    # align w/ BWA: -t for number of cores
                    bwa \
                      mem \
                      -R "@RG\tID:$(inputs.read_group_id)\tSM:$(inputs.read_group_sample_name)" \
                      -t $(inputs.cpu) \
                      $REF_DIR/genome.fa \
                      <(zcat $(inputs.fastq_input_read1.path)) <(zcat $(inputs.fastq_input_read2.path)) \
                      > $(inputs.output_base_name + ".aligned.sam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/snaptools:0.0.1
          - class: ResourceRequirement
            coresMin: $(inputs.cpu)
            ramMin: 3840.0
            outdirMin: '$((Math.ceil(3.25 * (function(size_of=0){inputs.fastq_input_read1.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.fastq_input_read2.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.tar_bwa_reference.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.fastq_input_read1.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.fastq_input_read2.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.tar_bwa_reference.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: SamToBam
    in:
      - id: sam_input
        source: BWAPairedEndAlignment/sam_aligned_output
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: SamToBam.docker_image
    out:
      - id: bam_output
      - id: monitoring_log
    run:
        id: SamToBam
        class: CommandLineTool
        inputs:
          - id: sam_input
            doc: the aligned sam produced by the aligner
            type: File
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using samtools to be used (default: quay.io/biocontainers/samtools:1.9--h10a08f8_12)'
            default: quay.io/biocontainers/samtools:1.9--h10a08f8_12
            type: string
        outputs:
          - id: bam_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".bam")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    # converst sam to bam
                    samtools view \
                      -bhS \
                      $(inputs.sam_input.path) \
                      -o $(inputs.output_base_name + ".bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/biocontainers/samtools:1.9--h10a08f8_12
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3840.0
            outdirMin: '$((Math.ceil(2 * (function(size_of=0){inputs.sam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.sam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: SortCoordinateOrder
    in:
      - id: bam_input
        source: SamToBam/bam_output
      - id: sort_order
        source: SortCoordinateOrder.sort_order
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: SortCoordinateOrder.docker_image
    out:
      - id: bam_sort_output
      - id: monitoring_log
    run:
        id: SortSam
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam to be sorted by picard tools
            type: File
          - id: sort_order
            doc: 'the desired way for the bam to be sorted (default: coordinate)'
            default: coordinate
            type: string
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using picard to be used (default: quay.io/biocontainers/picard:2.18.3--py27_0)'
            default: quay.io/biocontainers/picard:2.18.3--py27_0
            type: string
        outputs:
          - id: bam_sort_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".sorted." + inputs.sort_order +
                    ".bam")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    java -Xmx3250m -jar /usr/local/share/picard-2.18.3-0/picard.jar SortSam \
                      INPUT=$(inputs.bam_input.path) \
                      SORT_ORDER=$(inputs.sort_order) \
                      MAX_RECORDS_IN_RAM=300000 \
                      OUTPUT=$(inputs.output_base_name + ".sorted." + inputs.sort_order + ".bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/biocontainers/picard:2.18.3--py27_0
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3750.0
            outdirMin: '$((Math.ceil(3.25 * (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: FilterMinMapQuality
    in:
      - id: bam_input
        source: SortCoordinateOrder/bam_sort_output
      - id: min_map_quality
        source: min_map_quality
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: FilterMinMapQuality.docker_image
    out:
      - id: bam_filter_mapq_output
      - id: monitoring_log
    run:
        id: FilterMinMapQuality
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam to passed into samtools tools
            type: File
          - id: min_map_quality
            doc: the minimum mapping quality to be filtered by samtools view and snap-pre
                (snaptools task)
            type: int
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using samtools to be used (default: quay.io/biocontainers/samtools:1.9--h10a08f8_12)'
            default: quay.io/biocontainers/samtools:1.9--h10a08f8_12
            type: string
        outputs:
          - id: bam_filter_mapq_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".filtered.min_map_quality.bam")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    # filter for a map quality
                    # -b output is bam, -h include header, -q reads with mapping quality >=
                    samtools view \
                      -bh \
                      -q$(inputs.min_map_quality) \
                      $(inputs.bam_input.path) \
                      > $(inputs.output_base_name + ".filtered.min_map_quality.bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/biocontainers/samtools:1.9--h10a08f8_12
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3840.0
            outdirMin: '$((Math.ceil(2 * (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: FilterMaxFragmentLength
    in:
      - id: bam_input
        source: FilterMinMapQuality/bam_filter_mapq_output
      - id: max_fragment_length
        source: max_fragment_length
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: FilterMaxFragmentLength.docker_image
    out:
      - id: bam_filter_fragment_length_output
      - id: monitoring_log
    run:
        id: FilterMaxFragmentLength
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam to passed into gatk tools
            type: File
          - id: max_fragment_length
            doc: the maximum fragment length for filtering out reads by gatk (snaptools
                task)
            type: int
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using gatk to be used (default: broadinstitute/gatk:4.1.2.0)'
            default: broadinstitute/gatk:4.1.2.0
            type: string
        outputs:
          - id: bam_filter_fragment_length_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".filtered.max_fragment_length.bam")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    gatk --java-options "-Xms2750m -Xmx3250m" \
                      PrintReads \
                      --input=$(inputs.bam_input.path) \
                      --read-filter FragmentLengthReadFilter --max-fragment-length $(inputs.max_fragment_length) \
                      --output=$(inputs.output_base_name + ".filtered.max_fragment_length.bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: broadinstitute/gatk:4.1.2.0
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3750.0
            outdirMin: '$((Math.ceil(2 * (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: FilterMitochondrialReads
    in:
      - id: bam_input
        source: FilterMaxFragmentLength/bam_filter_fragment_length_output
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: FilterMitochondrialReads.docker_image
    out:
      - id: bam_no_chrM_reads_output
      - id: bam_chrM_reads_output
      - id: monitoring_log
    run:
        id: FilterMitochondrialReads
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam to passed into samtools tools
            type: File
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using samtools to be used (default: quay.io/biocontainers/samtools:1.9--h10a08f8_12)'
            default: quay.io/biocontainers/samtools:1.9--h10a08f8_12
            type: string
        outputs:
          - id: bam_no_chrM_reads_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".filtered.no_mitochondrial_reads.bam")
          - id: bam_chrM_reads_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".filtered.mitochondrial_reads.bam")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    # create bam index to filter by chromosome
                    mkdir bam_input
                    ln -s $(inputs.bam_input.path) bam_input/$(inputs.bam_input.basename)
                    samtools index -b bam_input/$(inputs.bam_input.basename)

                    #get list of chromosomes from bam header (ignoring chrM chromosome)
                    declare -r LIST_CHRS=`samtools view -H bam_input/$(inputs.bam_input.basename) \
                      | grep chr \
                      | cut -f2 \
                      | sed 's/SN://g' \
                      | grep -v 'chrM\|_'`

                    # get bam w/o chrM using the list
                    samtools view \
                      -bh \
                      -f 0x2 \
                      bam_input/$(inputs.bam_input.basename) \
                      `echo $LIST_CHRS` \
                      -o $(inputs.output_base_name + ".filtered.no_mitochondrial_reads.bam")

                    #get bam with only chrM
                    samtools view  \
                      -bh \
                      bam_input/$(inputs.bam_input.basename) \
                      chrM \
                      -o $(inputs.output_base_name + ".filtered.mitochondrial_reads.bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/biocontainers/samtools:1.9--h10a08f8_12
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3840.0
            outdirMin: '$((Math.ceil(2 * (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: MakeCompliantChrMBAM
    in:
      - id: bam_input
        source: FilterMitochondrialReads/bam_chrM_reads_output
      - id: output_base_name
        source: output_base_name
        valueFrom: $(self + ".chrM_reads")
      - id: docker_image
        source: MakeCompliantChrMBAM.docker_image
    out:
      - id: compliant_bam_output
    run:
        id: MakeCompliantBAM
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam with barcodes in the read ids that need to be converted to
                barcodes in bam tags
            type: File
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using the python script to convert the bam barcodes/read
                ids (default: quay.io/humancellatlas/snaptools:0.0.1)'
            default: quay.io/humancellatlas/snaptools:0.0.1
            type: string
        outputs:
          - id: compliant_bam_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".compliant.bam")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    makeCompliantBAM.py \
                      --input-bam $(inputs.bam_input.path) \
                      --output-bam $(inputs.output_base_name + ".compliant.bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/snaptools:0.0.1
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 4096.0
            outdirMin: '$((Math.ceil(2.5 * (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: SortQueryName
    in:
      - id: bam_input
        source: FilterMitochondrialReads/bam_no_chrM_reads_output
      - id: sort_order
        default: queryname
      - id: output_base_name
        source: output_base_name
      - id: docker_image
        source: SortQueryName.docker_image
    out:
      - id: bam_sort_output
      - id: monitoring_log
    run:
        id: SortSam
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam to be sorted by picard tools
            type: File
          - id: sort_order
            doc: 'the desired way for the bam to be sorted (default: coordinate)'
            default: coordinate
            type: string
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using picard to be used (default: quay.io/biocontainers/picard:2.18.3--py27_0)'
            default: quay.io/biocontainers/picard:2.18.3--py27_0
            type: string
        outputs:
          - id: bam_sort_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".sorted." + inputs.sort_order +
                    ".bam")
          - id: monitoring_log
            type:
              - File
              - 'null'
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    java -Xmx3250m -jar /usr/local/share/picard-2.18.3-0/picard.jar SortSam \
                      INPUT=$(inputs.bam_input.path) \
                      SORT_ORDER=$(inputs.sort_order) \
                      MAX_RECORDS_IN_RAM=300000 \
                      OUTPUT=$(inputs.output_base_name + ".sorted." + inputs.sort_order + ".bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/biocontainers/picard:2.18.3--py27_0
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3750.0
            outdirMin: '$((Math.ceil(3.25 * (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: MakeCompliantFilteredAndSortedBAM
    in:
      - id: bam_input
        source: SortQueryName/bam_sort_output
      - id: output_base_name
        source: output_base_name
        valueFrom: $(self + ".filtered_and_sorted")
      - id: docker_image
        source: MakeCompliantFilteredAndSortedBAM.docker_image
    out:
      - id: compliant_bam_output
    run:
        id: MakeCompliantBAM
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam with barcodes in the read ids that need to be converted to
                barcodes in bam tags
            type: File
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: docker_image
            doc: 'the docker image using the python script to convert the bam barcodes/read
                ids (default: quay.io/humancellatlas/snaptools:0.0.1)'
            default: quay.io/humancellatlas/snaptools:0.0.1
            type: string
        outputs:
          - id: compliant_bam_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".compliant.bam")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    makeCompliantBAM.py \
                      --input-bam $(inputs.bam_input.path) \
                      --output-bam $(inputs.output_base_name + ".compliant.bam")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/snaptools:0.0.1
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 4096.0
            outdirMin: '$((Math.ceil(2.5 * (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.bam_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: SnapPre
    in:
      - id: bam_input
        source: SortQueryName/bam_sort_output
      - id: output_base_name
        source: output_base_name
      - id: genome_name
        source: genome_name
      - id: max_fragment_length
        source: max_fragment_length
      - id: genome_size_file
        source: genome_size_file
      - id: docker_image
        source: SnapPre.docker_image
    out:
      - id: snap_file_output
      - id: snap_qc_output
    run:
        id: SnapPre
        class: CommandLineTool
        inputs:
          - id: bam_input
            doc: the bam to passed into snaptools tools
            type: File
          - id: output_base_name
            doc: base name to be used for the output of the task
            type: string
          - id: genome_name
            doc: the name of the genome being analyzed
            type: string
          - id: max_fragment_length
            doc: the maximum fragment length for filtering out reads by snap-pre (snaptools
                task)
            type: int
          - id: genome_size_file
            doc: 'size for the chromoomes for the genome; ex: mm10.chrom.size'
            type: File
          - id: docker_image
            doc: 'the docker image using snaptools to be used (default: quay.io/humancellatlas/snaptools:0.0.1)'
            default: quay.io/humancellatlas/snaptools:0.0.1
            type: string
        outputs:
          - id: snap_file_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".snap")
          - id: snap_qc_output
            type: File
            outputBinding:
                glob: $(inputs.output_base_name + ".snap" + ".qc")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    # Does the main counting
                    snaptools snap-pre \
                      --input-file=$(inputs.bam_input.path) \
                      --output-snap=$(inputs.output_base_name + ".snap") \
                      --genome-name=$(inputs.genome_name) \
                      --genome-size=$(inputs.genome_size_file.path) \
                      --min-mapq=0 \
                      --min-flen=0 \
                      --max-flen=$(inputs.max_fragment_length) \
                      --keep-chrm=TRUE \
                      --keep-single=TRUE \
                      --keep-secondary=False \
                      --overwrite=True \
                      --max-num=1000000 \
                      --min-cov=100 \
                      --verbose=True
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/snaptools:0.0.1
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 16384.0
            outdirMin: 153600
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: SnapCellByBin
    in:
      - id: snap_input
        source: SnapPre/snap_file_output
      - id: bin_size_list
        source: bin_size_list
      - id: snap_output_name
        source: SnapCellByBin.snap_output_name
      - id: docker_image
        source: SnapCellByBin.docker_image
    out:
      - id: snap_output
    run:
        id: SnapCellByBin
        class: CommandLineTool
        inputs:
          - id: snap_input
            doc: the bam to passed into snaptools tools
            type: File
          - id: bin_size_list
            doc: space separated list of bins to generate
            type: string
          - id: snap_output_name
            doc: output.snap
            default: output.snap
            type: string
          - id: docker_image
            doc: 'the docker image to be used (default: quay.io/humancellatlas/snaptools:0.0.1)'
            default: quay.io/humancellatlas/snaptools:0.0.1
            type: string
        outputs:
          - id: snap_output
            type: File
            outputBinding:
                glob: $(inputs.snap_output_name)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail

                    cp $(inputs.snap_input.path) $(inputs.snap_output_name)

                    # This is mutating the file in-place
                    snaptools snap-add-bmat  \
                      --snap-file=$(inputs.snap_output_name)  \
                      --bin-size-list $(inputs.bin_size_list)  \
                      --verbose=True
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/snaptools:0.0.1
          - class: ResourceRequirement
            coresMin: $(1)
            ramMin: 16384.0
            outdirMin: 153600
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BreakoutSnap
    in:
      - id: snap_input
        source: SnapCellByBin/snap_output
      - id: docker_image
        source: BreakoutSnap.docker_image
      - id: bin_size_list
        source: bin_size_list
    out:
      - id: barcodes
      - id: fragments
      - id: binCoordinates
      - id: binCounts
      - id: barcodesSection
    run:
        id: BreakoutSnap
        class: CommandLineTool
        inputs:
          - id: snap_input
            type: File
          - id: docker_image
            default: quay.io/humancellatlas/snap-breakout:0.0.1
            type: string
          - id: bin_size_list
            type: string
        outputs:
          - id: barcodes
            type: File
            outputBinding:
                glob: output/barcodes.csv
          - id: fragments
            type: File
            outputBinding:
                glob: output/fragments.csv
          - id: binCoordinates
            type: File
            outputBinding:
                glob: $('output/binCoordinates_' + inputs.bin_size_list + '.csv')
          - id: binCounts
            type: File
            outputBinding:
                glob: $('output/binCounts_' + inputs.bin_size_list + '.csv')
          - id: barcodesSection
            type: File
            outputBinding:
                glob: output/barcodesSection.csv
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -euo pipefail
                    mkdir output
                    breakoutSnap.py --input $(inputs.snap_input.path) \
                        --output-prefix output/
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/snap-breakout:0.0.1
          - class: ResourceRequirement
            coresMin: $(1)
            ramMin: 15258.7890625
            outdirMin: '$((Math.ceil(10 * (function(size_of=0){inputs.snap_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3 < 1 ? 1 : (function(size_of=0){inputs.snap_input.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1024^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
outputs:
  - id: ATAC.bam_chrM_reads_compliant_output
    outputSource: MakeCompliantChrMBAM/compliant_bam_output
    type: File
  - id: ATAC.bam_filtered_and_sorted_compliant_output
    outputSource: MakeCompliantFilteredAndSortedBAM/compliant_bam_output
    type: File
  - id: ATAC.snap_qc_output
    outputSource: SnapPre/snap_qc_output
    type: File
  - id: ATAC.snap_output
    outputSource: SnapCellByBin/snap_output
    type: File
