cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: CreateSequenceGroupingTSV
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                python3 <<CODE
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
        dockerPull: us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian
      - class: ResourceRequirement
        ramMin: 2048.0
        outdirMin: 1024
    inputs:
      - id: ref_dict
        type: File
      - id: preemptible_tries
        type: int
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
    id: ScatterIntervalList
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                mkdir out
                java -Xms1000m -Xmx1500m -jar /usr/gitc/picard.jar \
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
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z
      - class: ResourceRequirement
        ramMin: 2000.0
        outdirMin: 1024
    inputs:
      - id: interval_list
        type: File
      - id: scatter_count
        type: int
      - id: break_bands_at_multiples_of
        type: int
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
            entry: |2

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
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616
      - class: ResourceRequirement
        coresMin: 1
        ramMin: 3072.0
        outdirMin: $((Math.ceil(2 * (function(size_of=0){inputs.input_bam.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.path.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: output_basename
        type: string
      - id: preemptible_tries
        type: int
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
    id: ConvertToBam
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                set -e
                set -o pipefail

                samtools view -b -o $(inputs.output_basename).bam -T $(inputs.ref_fasta.path) $(inputs.input_cram.path)

                samtools index $(inputs.output_basename).bam
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616
      - class: ResourceRequirement
        coresMin: 1
        ramMin: 3072.0
        outdirMin: 204800
    inputs:
      - id: input_cram
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: output_basename
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: output_bam
        type: File
        outputBinding:
            glob: $(inputs.output_basename + ".bam")
      - id: output_bam_index
        type: File
        outputBinding:
            glob: $(inputs.output_basename + ".bam.bai")
  - cwlVersion: v1.2
    id: SumFloats
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                python3 -c 'print($(inputs.sizes.join("+")))'
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian
      - class: ResourceRequirement
        outdirMin: 1024
    inputs:
      - id: sizes
        type:
            items: float
            type: array
      - id: preemptible_tries
        type: int
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
  - cwlVersion: v1.2
    id: ErrorWithMessage
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |2

                >&2 echo "Error: $(inputs.message)"
                exit 1
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: ubuntu:20.04
      - class: ResourceRequirement
        outdirMin: 1024
    inputs:
      - id: message
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs: []
