cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: GetBwaVersion
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
                # the sed may also fail with that error and that is something we actually want to fail on.
                /usr/gitc/bwa 2>&1 | \
                grep -e '^Version' | \
                sed 's/Version: //'
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        ramMin: 953.67431640625
        outdirMin: 1024
    inputs: []
    baseCommand:
      - bash
      - script.bash
    stdout: _stdout
    outputs:
      - id: version
        type: string
        outputBinding:
            loadContents: true
            glob: _stdout
            outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
  - cwlVersion: v1.2
    id: SamToFastqAndBwaMemAndMba
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -o pipefail
                set -e

                # set the bash variable needed for the command-line
                bash_ref_fasta=$(inputs.ref_fasta.path)
                # if ref_alt has data in it,
                if [ -s $(inputs.ref_alt.path) ]; then
                  java -Xms5000m -jar /usr/gitc/picard.jar \
                    SamToFastq \
                    INPUT=$(inputs.input_bam.path) \
                    FASTQ=/dev/stdout \
                    INTERLEAVE=true \
                    NON_PF=true | \
                  /usr/gitc/$(inputs.bwa_commandline) /dev/stdin - 2> >(tee $(inputs.output_bam_basename).bwa.stderr.log >&2) | \
                  java -Dsamjdk.compression_level=$(inputs.compression_level) -Xms3000m -jar /usr/gitc/picard.jar \
                    MergeBamAlignment \
                    VALIDATION_STRINGENCY=SILENT \
                    EXPECTED_ORIENTATIONS=FR \
                    ATTRIBUTES_TO_RETAIN=X0 \
                    ATTRIBUTES_TO_REMOVE=NM \
                    ATTRIBUTES_TO_REMOVE=MD \
                    ALIGNED_BAM=/dev/stdin \
                    UNMAPPED_BAM=$(inputs.input_bam.path) \
                    OUTPUT=$(inputs.output_bam_basename).bam \
                    REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
                    PAIRED_RUN=true \
                    SORT_ORDER="unsorted" \
                    IS_BISULFITE_SEQUENCE=false \
                    ALIGNED_READS_ONLY=false \
                    CLIP_ADAPTERS=false \
                    MAX_RECORDS_IN_RAM=2000000 \
                    ADD_MATE_CIGAR=true \
                    MAX_INSERTIONS_OR_DELETIONS=-1 \
                    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
                    PROGRAM_RECORD_ID="bwamem" \
                    PROGRAM_GROUP_VERSION="$(inputs.bwa_version)" \
                    PROGRAM_GROUP_COMMAND_LINE="$(inputs.bwa_commandline)" \
                    PROGRAM_GROUP_NAME="bwamem" \
                    UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
                    ALIGNER_PROPER_PAIR_FLAGS=true \
                    UNMAP_CONTAMINANT_READS=true \
                    ADD_PG_TAG_TO_READS=false

                  grep -m1 "read .* ALT contigs" $(inputs.output_bam_basename).bwa.stderr.log | \
                  grep -v "read 0 ALT contigs"

                # else ref_alt is empty or could not be found
                else
                  exit 1;
                fi
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        coresMin: 16
        ramMin: 13351.4404296875
        outdirMin: $(parseFloat(inputs.disk_size.replace("\\..*", "") ) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: output_bam
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam")
      - id: bwa_stderr_log
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bwa.stderr.log")
  - cwlVersion: v1.2
    id: SortSam
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4+

                java -Dsamjdk.compression_level=$(inputs.compression_level) -Xms4000m -jar /usr/gitc/picard.jar \
                  SortSam \
                  INPUT=$(inputs.input_bam.path) \
                  OUTPUT=$(inputs.output_bam_basename).bam \
                  SORT_ORDER="coordinate" \
                  CREATE_INDEX=true \
                  CREATE_MD5_FILE=true \
                  MAX_RECORDS_IN_RAM=300000

      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        coresMin: 1
        ramMin: 4768.37158203125
        outdirMin: $(parseFloat(inputs.disk_size.replace("\\..*", "") ) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: output_bam
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam")
      - id: output_bam_index
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bai")
      - id: output_bam_md5
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam.md5")
  - cwlVersion: v1.2
    id: MarkDuplicates
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Dsamjdk.compression_level=$(inputs.compression_level) -Xms4000m -jar /usr/gitc/picard.jar \
                  MarkDuplicates \
                  INPUT=$(inputs.input_bams.map(function(el) {return el.path}).join(" INPUT=")) \
                  OUTPUT=$(inputs.output_bam_basename).bam \
                  METRICS_FILE=$(inputs.metrics_filename) \
                  VALIDATION_STRINGENCY=SILENT \
                  $(inputs.read_name_regex === null ? "" : "READ_NAME_REGEX=" + inputs.read_name_regex) \
                  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
                  ASSUME_SORT_ORDER="queryname" \
                  CLEAR_DT="false" \
                  ADD_PG_TAG_TO_READS=false
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        ramMin: 6675.72021484375
        outdirMin: $(parseFloat(inputs.disk_size.replace("\\..*", "") ) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: output_bam
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam")
      - id: duplicate_metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_filename)
  - cwlVersion: v1.2
    id: CreateSequenceGroupingTSV
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                python <<CODE
                with open("$(inputs.ref_dict.path)", "r") as ref_dict_file:
                    sequence_tuple_list = []
                    longest_sequence = 0
                    for line in ref_dict_file:
                        if line.startswith("@SQ"):
                            line_split = line.split("\t")
                            # (Sequence_Name, Sequence_Length)
                            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
                    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
                # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
                # the last element after a :, so we add this as a sacrificial element.
                hg38_protection_tag = ":1+"
                # initialize the tsv string with the first sequence
                tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
                temp_size = sequence_tuple_list[0][1]
                for sequence_tuple in sequence_tuple_list[1:]:
                    if temp_size + sequence_tuple[1] <= longest_sequence:
                        temp_size += sequence_tuple[1]
                        tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
                    else:
                        tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
                        temp_size = sequence_tuple[1]
                # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
                with open("sequence_grouping.txt","w") as tsv_file:
                  tsv_file.write(tsv_string)
                  tsv_file.close()

                tsv_string += '\n' + "unmapped"

                with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
                  tsv_file_with_unmapped.write(tsv_string)
                  tsv_file_with_unmapped.close()
                CODE
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: python:2.7
      - class: ResourceRequirement
        ramMin: 1907.3486328125
        outdirMin: 1024
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: sequence_grouping
        type:
            items:
                items: string
                type: array
            type: array
        outputBinding:
            loadContents: true
            glob: sequence_grouping.txt
            outputEval: |-
                ${
                  var result = Array();
                  var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
                  // ^ remove any trailing newline to prevent a null being returned
                  contents.split(/\r\n|\r|\n/).forEach(function(line) {
                    result.push(line.split('\t'));
                  });
                  return result;
                }
      - id: sequence_grouping_with_unmapped
        type:
            items:
                items: string
                type: array
            type: array
        outputBinding:
            loadContents: true
            glob: sequence_grouping_with_unmapped.txt
            outputEval: |-
                ${
                  var result = Array();
                  var contents = self[0].contents.replace(/\r\n$/, "").replace(/\n$/, "").replace(/\r$/, "");
                  // ^ remove any trailing newline to prevent a null being returned
                  contents.split(/\r\n|\r|\n/).forEach(function(line) {
                    result.push(line.split('\t'));
                  });
                  return result;
                }
  - cwlVersion: v1.2
    id: BaseRecalibrator
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                /usr/gitc/gatk4/gatk-launch --javaOptions "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
                  -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
                  -Xloggc:gc_log.log -Xms4000m" \
                  BaseRecalibrator \
                  -R $(inputs.ref_fasta.path) \
                  -I $(inputs.input_bam) \
                  --useOriginalQualities \
                  -O $(inputs.recalibration_report_filename) \
                  -knownSites $(inputs.dbSNP_vcf.path) \
                  -knownSites $(inputs.known_indels_sites_VCFs.map(function(el) {return el.path}).join(" -knownSites ")) \
                  -L $(inputs.sequence_group_interval.join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        ramMin: 5722.0458984375
        outdirMin: $(parseFloat(inputs.disk_size.replace("\\..*", "") ) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: recalibration_report
        type: File
        outputBinding:
            glob: $(inputs.recalibration_report_filename)
  - cwlVersion: v1.2
    id: ApplyBQSR
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                /usr/gitc/gatk4/gatk-launch --javaOptions "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
                  -XX:+PrintGCDetails -Xloggc:gc_log.log \
                  -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=$(inputs.compression_level) -Xms3000m" \
                  ApplyBQSR \
                  --createOutputBamMD5 \
                  --addOutputSAMProgramRecord \
                  -R $(inputs.ref_fasta.path) \
                  -I $(inputs.input_bam) \
                  --useOriginalQualities \
                  -O $(inputs.output_bam_basename).bam \
                  -bqsr $(inputs.recalibration_report.path) \
                  -SQQ 10 -SQQ 20 -SQQ 30 \
                  -L $(inputs.sequence_group_interval.join(" -L "))
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        ramMin: 3337.860107421875
        outdirMin: $(parseFloat(inputs.disk_size.replace("\\..*", "") ) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: recalibrated_bam
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam")
      - id: recalibrated_bam_checksum
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam.md5")
  - cwlVersion: v1.2
    id: GatherBqsrReports
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                /usr/gitc/gatk4/gatk-launch --javaOptions "-Xms3000m" \
                  GatherBQSRReports \
                  -I $(inputs.input_bqsr_reports.map(function(el) {return el.path}).join(" -I ")) \
                  -O $(inputs.output_report_filename)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        ramMin: 3337.860107421875
        outdirMin: $((inputs.disk_size) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: output_bqsr_report
        type: File
        outputBinding:
            glob: $(inputs.output_report_filename)
  - cwlVersion: v1.2
    id: GatherBamFiles
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Dsamjdk.compression_level=$(inputs.compression_level) -Xms2000m -jar /usr/gitc/picard.jar \
                  GatherBamFiles \
                  INPUT=$(inputs.input_bams.map(function(el) {return el.path}).join(" INPUT=")) \
                  OUTPUT=$(inputs.output_bam_basename).bam \
                  CREATE_INDEX=true \
                  CREATE_MD5_FILE=true
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        ramMin: 2861.02294921875
        outdirMin: $(parseFloat(inputs.disk_size.replace("\\..*", "") ) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: output_bam
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam")
      - id: output_bam_index
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bai")
      - id: output_bam_md5
        type: File
        outputBinding:
            glob: $(inputs.output_bam_basename + ".bam.md5")
  - cwlVersion: v1.2
    id: ScatterIntervalList
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir out
                java -Xms1g -jar /usr/gitc/picard.jar \
                  IntervalListTools \
                  SCATTER_COUNT=$(inputs.scatter_count) \
                  SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
                  UNIQUE=true \
                  SORT=true \
                  BREAK_BANDS_AT_MULTIPLES_OF=$(inputs.break_bands_at_multiples_of) \
                  INPUT=$(inputs.interval_list.path) \
                  OUTPUT=out

                python3 <<CODE
                import glob, os
                # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
                intervals = sorted(glob.glob("out/*/*.interval_list"))
                for i, interval in enumerate(intervals):
                  (directory, filename) = os.path.split(interval)
                  newName = os.path.join(directory, str(i + 1) + filename)
                  os.rename(interval, newName)
                print(len(intervals))
                CODE
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        ramMin: 1907.3486328125
        outdirMin: 1024
    inputs: []
    baseCommand:
      - bash
      - script.bash
    stdout: _stdout
    outputs:
      - id: out
        type:
            items: File
            type: array
        outputBinding:
            glob: $("out/*/*.interval_list")
      - id: interval_count
        type: int
        outputBinding:
            loadContents: true
            glob: _stdout
            outputEval: $(parseInt(self[0].contents))
  - cwlVersion: v1.2
    id: ConvertToCram
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                set -o pipefail

                samtools view -C -T $(inputs.ref_fasta.path) $(inputs.input_bam.path) | \
                tee $(inputs.output_basename).cram | \
                md5sum | awk '{print $1}' > $(inputs.output_basename).cram.md5

                # Create REF_CACHE. Used when indexing a CRAM
                seq_cache_populate.pl -root ./ref/cache $(inputs.ref_fasta.path)
                export REF_PATH=:
                export REF_CACHE=./ref/cache/%2s/%2s/%s

                samtools index $(inputs.output_basename).cram
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: ResourceRequirement
        coresMin: 1
        ramMin: 2861.02294921875
        outdirMin: $(parseFloat(inputs.disk_size.replace("\\..*", "") ) * 1024)
    inputs: []
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: output_cram
        type: File
        outputBinding:
            glob: $(inputs.output_basename + ".cram")
      - id: output_cram_index
        type: File
        outputBinding:
            glob: $(inputs.output_basename + ".cram.crai")
      - id: output_cram_md5
        type: File
        outputBinding:
            glob: $(inputs.output_basename + ".cram.md5")
  - cwlVersion: v1.2
    id: SumFloats
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                python -c "print $()"
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: python:2.7
      - class: ResourceRequirement
        outdirMin: 1024
    inputs: []
    baseCommand:
      - bash
      - script.bash
    stdout: _stdout
    outputs:
      - id: total_size
        type: float
        outputBinding:
            loadContents: true
            glob: _stdout
            outputEval: $(parseFloat(self[0].contents))
