cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: CollectQualityYieldMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms2000m -Xmx3000m -jar /usr/picard/picard.jar \
                  CollectQualityYieldMetrics \
                  INPUT=$(inputs.input_bam.path) \
                  OQ=true \
                  OUTPUT=$(inputs.metrics_filename)
                sed -i -e 1,5d $(inputs.metrics_filename)  # for reproducibility
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 3500.0
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: metrics_filename
        type: string
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: quality_yield_metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_filename)
  - cwlVersion: v1.2
    id: CollectUnsortedReadgroupBamQualityMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
                  CollectMultipleMetrics \
                  INPUT=$(inputs.input_bam.path) \
                  OUTPUT=$(inputs.output_bam_prefix) \
                  ASSUME_SORTED=true \
                  PROGRAM=null \
                  PROGRAM=CollectBaseDistributionByCycle \
                  PROGRAM=CollectInsertSizeMetrics \
                  PROGRAM=MeanQualityByCycle \
                  PROGRAM=QualityScoreDistribution \
                  METRIC_ACCUMULATION_LEVEL=null \
                  METRIC_ACCUMULATION_LEVEL=ALL_READS

                touch $(inputs.output_bam_prefix).insert_size_metrics
                touch $(inputs.output_bam_prefix).insert_size_histogram.pdf
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 7000.0
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: output_bam_prefix
        type: string
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: base_distribution_by_cycle_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".base_distribution_by_cycle.pdf")
      - id: base_distribution_by_cycle_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".base_distribution_by_cycle_metrics")
      - id: insert_size_histogram_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".insert_size_histogram.pdf")
      - id: insert_size_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".insert_size_metrics")
      - id: quality_by_cycle_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".quality_by_cycle.pdf")
      - id: quality_by_cycle_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".quality_by_cycle_metrics")
      - id: quality_distribution_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".quality_distribution.pdf")
      - id: quality_distribution_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".quality_distribution_metrics")
  - cwlVersion: v1.2
    id: CollectReadgroupBamQualityMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                # These are optionally generated, but need to exist for Cromwell's sake
                touch $(inputs.output_bam_prefix).gc_bias.detail_metrics \
                  $(inputs.output_bam_prefix).gc_bias.pdf \
                  $(inputs.output_bam_prefix).gc_bias.summary_metrics

                java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
                  CollectMultipleMetrics \
                  INPUT=$(inputs.input_bam.path) \
                  REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
                  OUTPUT=$(inputs.output_bam_prefix) \
                  ASSUME_SORTED=true \
                  PROGRAM=null \
                  PROGRAM=CollectAlignmentSummaryMetrics \
                  $(inputs.collect_gc_bias_metrics ? 'PROGRAM="CollectGcBiasMetrics"' : "") \
                  METRIC_ACCUMULATION_LEVEL=null \
                  METRIC_ACCUMULATION_LEVEL=READ_GROUP
                sed -i -e 1,5d "$(inputs.output_bam_prefix).alignment_summary_metrics"   # for reproducibility
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 7000.0
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_dict.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: input_bam_index
        type: File
      - id: output_bam_prefix
        type: string
      - id: ref_dict
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: collect_gc_bias_metrics
        default: true
        type: boolean
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: alignment_summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".alignment_summary_metrics")
      - id: gc_bias_detail_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".gc_bias.detail_metrics")
      - id: gc_bias_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".gc_bias.pdf")
      - id: gc_bias_summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".gc_bias.summary_metrics")
  - cwlVersion: v1.2
    id: CollectAggregationMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                # These are optionally generated, but need to exist for Cromwell's sake
                touch $(inputs.output_bam_prefix).gc_bias.detail_metrics \
                  $(inputs.output_bam_prefix).gc_bias.pdf \
                  $(inputs.output_bam_prefix).gc_bias.summary_metrics \
                  $(inputs.output_bam_prefix).insert_size_metrics \
                  $(inputs.output_bam_prefix).insert_size_histogram.pdf

                java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
                  CollectMultipleMetrics \
                  INPUT=$(inputs.input_bam.path) \
                  REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
                  OUTPUT=$(inputs.output_bam_prefix) \
                  ASSUME_SORTED=true \
                  PROGRAM=null \
                  PROGRAM=CollectAlignmentSummaryMetrics \
                  PROGRAM=CollectInsertSizeMetrics \
                  PROGRAM=CollectSequencingArtifactMetrics \
                  PROGRAM=QualityScoreDistribution \
                  $(inputs.collect_gc_bias_metrics ? 'PROGRAM="CollectGcBiasMetrics"' : "") \
                  METRIC_ACCUMULATION_LEVEL=null \
                  METRIC_ACCUMULATION_LEVEL=SAMPLE \
                  METRIC_ACCUMULATION_LEVEL=LIBRARY
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 7000.0
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_dict.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: input_bam_index
        type: File
      - id: output_bam_prefix
        type: string
      - id: ref_dict
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: collect_gc_bias_metrics
        default: true
        type: boolean
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: alignment_summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".alignment_summary_metrics")
      - id: bait_bias_detail_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".bait_bias_detail_metrics")
      - id: bait_bias_summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".bait_bias_summary_metrics")
      - id: gc_bias_detail_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".gc_bias.detail_metrics")
      - id: gc_bias_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".gc_bias.pdf")
      - id: gc_bias_summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".gc_bias.summary_metrics")
      - id: insert_size_histogram_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".insert_size_histogram.pdf")
      - id: insert_size_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".insert_size_metrics")
      - id: pre_adapter_detail_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".pre_adapter_detail_metrics")
      - id: pre_adapter_summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".pre_adapter_summary_metrics")
      - id: quality_distribution_pdf
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".quality_distribution.pdf")
      - id: quality_distribution_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".quality_distribution_metrics")
      - id: error_summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_bam_prefix + ".error_summary_metrics")
  - cwlVersion: v1.2
    id: ConvertSequencingArtifactToOxoG
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                input_base=\$(dirname $(inputs.pre_adapter_detail_metrics.path))/$(inputs.base_name)
                java -Xms$(Math.ceil(4000 * inputs.memory_multiplier)  - 1000)m -Xmx$(Math.ceil(4000 * inputs.memory_multiplier)  - 500)m \
                  -jar /usr/picard/picard.jar \
                  ConvertSequencingArtifactToOxoG \
                  --INPUT_BASE $input_base \
                  --OUTPUT_BASE $(inputs.base_name) \
                  --REFERENCE_SEQUENCE $(inputs.ref_fasta.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "MiB";
            var value = parseInt(`${Math.ceil(4000 * inputs.memory_multiplier) }`.match(/[0-9]+/g));
            var memory = "";
            if(unit==="KiB") memory = value/1024;
            else if(unit==="MiB") memory = value;
            else if(unit==="GiB") memory = value*1024;
            else if(unit==="TiB") memory = value*1024*1024;
            else if(unit==="B") memory = value/(1024*1024);
            else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
            else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
            else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
            else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: $((Math.ceil((function(size_of=0){inputs.pre_adapter_detail_metrics.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.bait_bias_detail_metrics.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_dict.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: pre_adapter_detail_metrics
        type: File
      - id: bait_bias_detail_metrics
        type: File
      - id: base_name
        type: string
      - id: ref_dict
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: preemptible_tries
        type: int
      - id: memory_multiplier
        default: 1
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: oxog_metrics
        type: File
        outputBinding:
            glob: $(inputs.base_name + ".oxog_metrics")
  - cwlVersion: v1.2
    id: CrossCheckFingerprints
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Dsamjdk.buffer_size=131072 \
                  -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx3000m \
                  -jar /usr/picard/picard.jar \
                  CrosscheckFingerprints \
                  OUTPUT=$(inputs.metrics_filename) \
                  HAPLOTYPE_MAP=$(inputs.haplotype_database_file.path) \
                  EXPECT_ALL_GROUPS_TO_MATCH=true \
                  INPUT=$(inputs.input_bams.map(function(el) {return el.path}).join(" INPUT=")) \
                  LOD_THRESHOLD=$(inputs.lod_threshold) \
                  CROSSCHECK_BY=$(inputs.cross_check_by)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 3500.0
        outdirMin: $((Math.ceil(inputs.total_input_size)  + 20) * 1024)
    inputs:
      - id: input_bams
        type:
            items: File
            type: array
      - id: input_bam_indexes
        type:
            items: File
            type: array
      - id: haplotype_database_file
        type: File
      - id: metrics_filename
        type: string
      - id: total_input_size
        type: float
      - id: preemptible_tries
        type: int
      - id: lod_threshold
        type: float
      - id: cross_check_by
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: cross_check_fingerprints_metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_filename)
  - cwlVersion: v1.2
    id: CheckFingerprint
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                java -Xms$(inputs.memory_size - 1000)m -Xmx$(inputs.memory_size - 500)m -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
                CheckFingerprint \
                  --INPUT $([inputs.input_vcf === null ? "" : inputs.input_vcf.path, inputs.input_bam === null ? "" : inputs.input_bam.path].find(function(element) { return element !== null }) ) \
                  $(inputs.input_vcf !== null ? inputs.input_sample_alias === null ? "" : "--OBSERVED_SAMPLE_ALIAS \"" + inputs.input_sample_alias + "\"" : "") \
                  --GENOTYPES $(inputs.genotypes.path) \
                  --EXPECTED_SAMPLE_ALIAS "$(inputs.expected_sample_alias)" \
                  $(inputs.input_bam !== null ? "--IGNORE_READ_GROUPS true" : "") \
                  --HAPLOTYPE_MAP $(inputs.haplotype_database_file.path) \
                  --GENOTYPE_LOD_THRESHOLD $(inputs.genotype_lod_threshold) \
                  --SUMMARY_OUTPUT $(inputs.output_basename + ".fingerprinting_summary_metrics") \
                  --DETAIL_OUTPUT $(inputs.output_basename + ".fingerprinting_detail_metrics") \
                  $(inputs.ref_fasta === null ? "" : "--REFERENCE_SEQUENCE " + inputs.ref_fasta.path)

                CONTENT_LINE=\$(cat $(inputs.output_basename + ".fingerprinting_summary_metrics") |
                grep -n "## METRICS CLASS\tpicard.analysis.FingerprintingSummaryMetrics" |
                cut -f1 -d:)
                CONTENT_LINE=\$(($CONTENT_LINE+2))
                sed '8q;d' $(inputs.output_basename + ".fingerprinting_summary_metrics") | cut -f5 > lod
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.26.4
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "MiB";
            var value = parseInt(`${inputs.memory_size}`.match(/[0-9]+/g));
            var memory = "";
            if(unit==="KiB") memory = value/1024;
            else if(unit==="MiB") memory = value;
            else if(unit==="GiB") memory = value*1024;
            else if(unit==="TiB") memory = value*1024*1024;
            else if(unit==="B") memory = value/(1024*1024);
            else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
            else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
            else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
            else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: '$((Math.ceil((function(size_of=0){inputs.input_bam === null ?
            "" : inputs.input_bam.path.forEach(function(element){ if (element) {size_of
            += element.size}})}) / 1024^3 + (function(size_of=0){inputs.input_vcf
            === null ? "" : inputs.input_vcf.path.forEach(function(element){ if (element)
            {size_of += element.size}})}) / 1024^3)  + 20) * 1024)'
    inputs:
      - id: input_bam
        type:
          - File
          - 'null'
      - id: input_bam_index
        type:
          - File
          - 'null'
      - id: input_vcf
        type:
          - File
          - 'null'
      - id: input_vcf_index
        type:
          - File
          - 'null'
      - id: input_sample_alias
        type:
          - string
          - 'null'
      - id: genotypes
        type: File
      - id: genotypes_index
        type:
          - File
          - 'null'
      - id: expected_sample_alias
        type: string
      - id: output_basename
        type: string
      - id: genotype_lod_threshold
        default: 5.0
        type: float
      - id: haplotype_database_file
        type: File
      - id: ref_fasta
        type:
          - File
          - 'null'
      - id: ref_fasta_index
        type:
          - File
          - 'null'
      - id: memory_size
        default: 2500
        type: int
      - id: preemptible_tries
        default: 3
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_basename + ".fingerprinting_summary_metrics")
      - id: detail_metrics
        type: File
        outputBinding:
            glob: $(inputs.output_basename + ".fingerprinting_detail_metrics")
      - id: lod
        type: float
        outputBinding:
            loadContents: true
            glob: lod
            outputEval: $(parseFloat(self[0].contents))
  - cwlVersion: v1.2
    id: CheckPreValidation
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4+

                set -o pipefail
                set -e

                grep -A 1 PERCENT_DUPLICATION $(inputs.duplication_metrics.path) > duplication.csv
                grep -A 3 PCT_CHIMERAS $(inputs.chimerism_metrics.path) | grep -v OF_PAIR > chimerism.csv

                python3 <<CODE

                import csv
                with open('duplication.csv') as dupfile:
                  reader = csv.DictReader(dupfile, delimiter='\t')
                  for row in reader:
                    with open("duplication_value.txt","w") as file:
                      file.write(row['PERCENT_DUPLICATION'])
                      file.close()

                with open('chimerism.csv') as chimfile:
                  reader = csv.DictReader(chimfile, delimiter='\t')
                  for row in reader:
                    with open("chimerism_value.txt","w") as file:
                      file.write(row['PCT_CHIMERAS'])
                      file.close()

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
      - id: duplication_metrics
        type: File
      - id: chimerism_metrics
        type: File
      - id: max_duplication_in_reasonable_sample
        type: float
      - id: max_chimerism_in_reasonable_sample
        type: float
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: duplication_rate
        type: float
        outputBinding:
            loadContents: true
            glob: duplication_value.txt
            outputEval: $(parseFloat(self[0].contents))
      - id: chimerism_rate
        type: float
        outputBinding:
            loadContents: true
            glob: chimerism_value.txt
            outputEval: $(parseFloat(self[0].contents))
      - id: is_outlier_data
        type: boolean
        outputBinding:
            outputEval: $("duplication_value.txt" > inputs.max_duplication_in_reasonable_sample
                || "chimerism_value.txt" > inputs.max_chimerism_in_reasonable_sample)
  - cwlVersion: v1.2
    id: ValidateSamFile
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms$(inputs.memory_size - 1000)m -Xmx$(inputs.memory_size - 500)m -jar /usr/picard/picard.jar \
                  ValidateSamFile \
                  INPUT=$(inputs.input_bam.path) \
                  OUTPUT=$(inputs.report_filename) \
                  REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
                  $(inputs.max_output === null ? "" : "MAX_OUTPUT=" + inputs.max_output) \
                  IGNORE=$(inputs.ignore === null ? "null" : inputs.ignore.join(" IGNORE=")) \
                  MODE=VERBOSE \
                  $(inputs.is_outlier_data === null ? "SKIP_MATE_VALIDATION=false" : inputs.is_outlier_data ? "SKIP_MATE_VALIDATION=true" : "SKIP_MATE_VALIDATION=false") \
                  IS_BISULFITE_SEQUENCED=false
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "MiB";
            var value = parseInt(`${inputs.memory_size}`.match(/[0-9]+/g));
            var memory = "";
            if(unit==="KiB") memory = value/1024;
            else if(unit==="MiB") memory = value;
            else if(unit==="GiB") memory = value*1024;
            else if(unit==="TiB") memory = value*1024*1024;
            else if(unit==="B") memory = value/(1024*1024);
            else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
            else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
            else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
            else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_dict.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + inputs.additional_disk)
            * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: input_bam_index
        type:
          - File
          - 'null'
      - id: report_filename
        type: string
      - id: ref_dict
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: max_output
        type:
          - int
          - 'null'
      - id: ignore
        type:
          - items: string
            type: array
          - 'null'
      - id: is_outlier_data
        type:
          - boolean
          - 'null'
      - id: preemptible_tries
        type: int
      - id: memory_multiplier
        default: 1
        type: int
      - id: additional_disk
        default: 20
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: report
        type: File
        outputBinding:
            glob: $(inputs.report_filename)
  - cwlVersion: v1.2
    id: CollectWgsMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
                  CollectWgsMetrics \
                  INPUT=$(inputs.input_bam.path) \
                  VALIDATION_STRINGENCY=SILENT \
                  REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
                  INCLUDE_BQ_HISTOGRAM=true \
                  INTERVALS=$(inputs.wgs_coverage_interval_list.path) \
                  OUTPUT=$(inputs.metrics_filename) \
                  USE_FAST_ALGORITHM=true \
                  READ_LENGTH=$(inputs.read_length)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 3000.0
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: input_bam_index
        type: File
      - id: metrics_filename
        type: string
      - id: wgs_coverage_interval_list
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: read_length
        type: int
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_filename)
  - cwlVersion: v1.2
    id: CollectRawWgsMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms$(inputs.memory_size - 1 * 1000)m -jar /usr/picard/picard.jar \
                  CollectRawWgsMetrics \
                  INPUT=$(inputs.input_bam.path) \
                  VALIDATION_STRINGENCY=SILENT \
                  REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
                  INCLUDE_BQ_HISTOGRAM=true \
                  INTERVALS=$(inputs.wgs_coverage_interval_list.path) \
                  OUTPUT=$(inputs.metrics_filename) \
                  USE_FAST_ALGORITHM=true \
                  READ_LENGTH=$(inputs.read_length)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "GiB";
            var value = parseInt(`${inputs.memory_size}`.match(/[0-9]+/g));
            var memory = "";
            if(unit==="KiB") memory = value/1024;
            else if(unit==="MiB") memory = value;
            else if(unit==="GiB") memory = value*1024;
            else if(unit==="TiB") memory = value*1024*1024;
            else if(unit==="B") memory = value/(1024*1024);
            else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
            else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
            else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
            else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + inputs.additional_disk)
            * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: input_bam_index
        type: File
      - id: metrics_filename
        type: string
      - id: wgs_coverage_interval_list
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: read_length
        type: int
      - id: preemptible_tries
        type: int
      - id: memory_multiplier
        default: 1
        type: int
      - id: additional_disk
        default: 20
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_filename)
  - cwlVersion: v1.2
    id: CollectHsMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms$(inputs.memory_size - 1000)m -Xmx$(inputs.memory_size - 500)m -jar /usr/picard/picard.jar \
                  CollectHsMetrics \
                  INPUT=$(inputs.input_bam.path) \
                  REFERENCE_SEQUENCE=$(inputs.ref_fasta.path) \
                  VALIDATION_STRINGENCY=SILENT \
                  TARGET_INTERVALS=$(inputs.target_interval_list.path) \
                  BAIT_INTERVALS=$(inputs.bait_interval_list.path) \
                  METRIC_ACCUMULATION_LEVEL=null \
                  METRIC_ACCUMULATION_LEVEL=SAMPLE \
                  METRIC_ACCUMULATION_LEVEL=LIBRARY \
                  OUTPUT=$(inputs.metrics_filename)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = "MiB";
            var value = parseInt(`${inputs.memory_size}`.match(/[0-9]+/g));
            var memory = "";
            if(unit==="KiB") memory = value/1024;
            else if(unit==="MiB") memory = value;
            else if(unit==="GiB") memory = value*1024;
            else if(unit==="TiB") memory = value*1024*1024;
            else if(unit==="B") memory = value/(1024*1024);
            else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
            else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
            else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
            else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
            else throw "Unknown units: " + unit;
            return parseInt(memory);
            }
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + inputs.additional_disk)
            * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: input_bam_index
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: metrics_filename
        type: string
      - id: target_interval_list
        type: File
      - id: bait_interval_list
        type: File
      - id: preemptible_tries
        type: int
      - id: memory_multiplier
        default: 1
        type: int
      - id: additional_disk
        default: 20
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_filename)
  - cwlVersion: v1.2
    id: CalculateReadGroupChecksum
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms1000m -Xmx3500m -jar /usr/picard/picard.jar \
                  CalculateReadGroupChecksum \
                  INPUT=$(inputs.input_bam.path) \
                  OUTPUT=$(inputs.read_group_md5_filename)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 4000.0
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_bam.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 40) * 1024)
    inputs:
      - id: input_bam
        type: File
      - id: input_bam_index
        type: File
      - id: read_group_md5_filename
        type: string
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: md5_file
        type: File
        outputBinding:
            glob: $(inputs.read_group_md5_filename)
  - cwlVersion: v1.2
    id: ValidateVCF
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                # Note that WGS needs a lot of memory to do the -L *.vcf if an interval file is not supplied
                gatk --java-options "-Xms6000m -Xmx6500m" \
                  ValidateVariants \
                  -V $(inputs.input_vcf.path) \
                  -R $(inputs.ref_fasta.path) \
                  -L $(inputs.calling_interval_list.path) \
                  $(inputs.is_gvcf ? "-gvcf" : "") \
                  --validation-type-to-exclude ALLELES \
                  $(inputs.dbsnp_vcf === null ? "" : "--dbsnp " + inputs.dbsnp_vcf.path) \
                  $(inputs.extra_args)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gatk/gatk:4.1.8.0
      - class: ResourceRequirement
        ramMin: 7000.0
        outdirMin: '$((Math.ceil((function(size_of=0){inputs.input_vcf.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.dbsnp_vcf
            === null ? "" : inputs.dbsnp_vcf.path.forEach(function(element){ if (element)
            {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_fasta_index.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.ref_dict.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)'
    inputs:
      - id: input_vcf
        type: File
      - id: input_vcf_index
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: ref_dict
        type: File
      - id: dbsnp_vcf
        type:
          - File
          - 'null'
      - id: dbsnp_vcf_index
        type:
          - File
          - 'null'
      - id: calling_interval_list
        type: File
      - id: calling_interval_list_index
        type:
          - File
          - 'null'
      - id: preemptible_tries
        default: 3
        type: int
      - id: is_gvcf
        default: true
        type: boolean
      - id: extra_args
        type:
          - string
          - 'null'
      - id: gatk_docker
        default: us.gcr.io/broad-gatk/gatk:4.1.8.0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs: []
  - cwlVersion: v1.2
    id: CollectVariantCallingMetrics
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
                  CollectVariantCallingMetrics \
                  INPUT=$(inputs.input_vcf.path) \
                  OUTPUT=$(inputs.metrics_basename) \
                  DBSNP=$(inputs.dbsnp_vcf.path) \
                  SEQUENCE_DICTIONARY=$(inputs.ref_dict.path) \
                  TARGET_INTERVALS=$(inputs.evaluation_interval_list.path) \
                  $(inputs.is_gvcf ? "GVCF_INPUT=true" : "")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8
      - class: ResourceRequirement
        ramMin: 3000.0
        outdirMin: $((Math.ceil((function(size_of=0){inputs.input_vcf.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3 + (function(size_of=0){inputs.dbsnp_vcf.forEach(function(element){
            if (element) {size_of += element.size}})}) / 1024^3)  + 20) * 1024)
    inputs:
      - id: input_vcf
        type: File
      - id: input_vcf_index
        type: File
      - id: metrics_basename
        type: string
      - id: dbsnp_vcf
        type: File
      - id: dbsnp_vcf_index
        type: File
      - id: ref_dict
        type: File
      - id: evaluation_interval_list
        type: File
      - id: is_gvcf
        default: true
        type: boolean
      - id: preemptible_tries
        type: int
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: summary_metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_basename + ".variant_calling_summary_metrics")
      - id: detail_metrics
        type: File
        outputBinding:
            glob: $(inputs.metrics_basename + ".variant_calling_detail_metrics")
