class: Workflow
id: mergeSvs
inputs:
  - id: estimate_sv_distance
    type: boolean
  - id: genome_build
    type: string
  - id: max_distance_to_merge
    type: int
  - id: minimum_sv_calls
    type: int
  - id: minimum_sv_size
    type: int
  - id: same_strand
    type: boolean
  - id: same_type
    type: boolean
  - id: snps_vcf
    type:
      - File
      - 'null'
  - id: sv_vcfs
    type:
        items: File
        type: array
  - id: blocklist_bedpe
    type:
      - File
      - 'null'
  - id: filterBlocklistSurvivor.slope
    default: 100
    type: int
  - id: bcftoolsMergeSvVcfs.force_merge
    default: true
    type: boolean
  - id: bcftoolsMergeSvVcfs.missing_ref
    default: false
    type: boolean
  - id: filterBlocklistBcftools.slope
    default: 100
    type: int
  - id: bcftoolsAnnotsvFilter.all_CDS
    default: false
    type: boolean
  - id: bcftoolsAnnotsvFilter.ignore_pass_filter
    default: false
    type: boolean
  - id: bcftoolsAnnotsvFilter.output_tsv_name
    default: filtered-bcftools-merged-AnnotSV.tsv
    type: string
outputs:
  - id: mergeSvs.bcftools_merged_sv_vcf
    outputSource: filterBlocklistBcftools/filtered_sv_vcf
    type: File
  - id: mergeSvs.bcftools_merged_annotated_tsv
    outputSource: bcftoolsAnnotateVariants/sv_variants_tsv
    type: File
  - id: mergeSvs.bcftools_merged_filtered_annotated_tsv
    outputSource: bcftoolsAnnotsvFilter/filtered_tsv
    type: File
  - id: mergeSvs.survivor_merged_sv_vcf
    outputSource: filterBlocklistSurvivor/filtered_sv_vcf
    type: File
  - id: mergeSvs.survivor_merged_annotated_tsv
    outputSource: survivorAnnotateVariants/sv_variants_tsv
    type: File
cwlVersion: v1.2
steps:
  - id: survivorMergeSvVcfs
    in:
      - id: vcfs
        source: sv_vcfs
      - id: max_distance_to_merge
        source: max_distance_to_merge
      - id: minimum_sv_calls
        source: minimum_sv_calls
      - id: same_type
        source: same_type
      - id: same_strand
        source: same_strand
      - id: estimate_sv_distance
        source: estimate_sv_distance
      - id: minimum_sv_size
        source: minimum_sv_size
      - id: cohort_name
        default: SURVIVOR-sv-merged/vcf
    out:
      - id: merged_vcf
    run:
        class: CommandLineTool
        id: survivor
        inputs:
          - id: vcfs
            type:
                items: File
                type: array
          - id: max_distance_to_merge
            type: int
          - id: minimum_sv_calls
            type: int
          - id: same_type
            type: boolean
          - id: same_strand
            type: boolean
          - id: estimate_sv_distance
            type: boolean
          - id: minimum_sv_size
            type: int
          - id: cohort_name
            default: SURVIVOR-sv-merged.vcf
            type: string
        outputs:
          - id: merged_vcf
            type: File
            outputBinding:
                glob: $(inputs.cohort_name)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    /bin/bash /usr/bin/survivor_merge_helper.sh \
                    $(inputs.vcfs.map(function(el) {return el.path}).join(",")) $(inputs.max_distance_to_merge) $(inputs.minimum_sv_calls) \
                    $(inputs.same_type ? "1" : "0") \
                    $(inputs.same_strand ? "1" : "0") \
                    $(inputs.estimate_sv_distance ? "1" : "0") \
                    $(inputs.minimum_sv_size) $(inputs.cohort_name)
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: mgibio/survivor-cwl:1.0.6.2
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 1907.3486328125
            outdirMin: 10240
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: filterBlocklistSurvivor
    in:
      - id: input_vcf
        source: survivorMergeSvVcfs/merged_vcf
      - id: blocklist_bedpe
        source: blocklist_bedpe
      - id: slope
        source: filterBlocklistSurvivor.slope
      - id: output_vcf_basename
        default: SURVIVOR-sv-merged
    out:
      - id: filtered_sv_vcf
      - id: filtered_sv_vcf_tbi
    run:
        class: CommandLineTool
        id: filterSvVcfBlocklistBedpe
        inputs:
          - id: input_vcf
            type: File
          - id: blocklist_bedpe
            type:
              - File
              - 'null'
          - id: slope
            default: 100
            type: int
          - id: output_vcf_basename
            default: blocklist_filtered
            type: string
        outputs:
          - id: filtered_sv_vcf
            type: File
            outputBinding:
                glob: $(inputs.output_vcf_basename + ".vcf.gz")
          - id: filtered_sv_vcf_tbi
            type: File
            outputBinding:
                glob: $(inputs.output_vcf_basename + ".vcf.gz.tbi")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eou pipefail
                    set -o errexit

                    INPUT_VCF="$(inputs.input_vcf.path)"
                    SLOPE="$(inputs.slope)"
                    OUT_BASE="$(inputs.output_vcf_basename)"
                    if [[ $(inputs.blocklist_bedpe) = true ]];then # blocklist_bedpe is passed.
                        BL_BEDPE="$(inputs.blocklist_bedpe === null ? "" : inputs.blocklist_bedpe.path)"

                        #CNVkit outputs invalid format like CIPOS=.,894;CIEND=.,894, which can cause svtools vcftobedpe fail
                        if [[ "$INPUT_VCF" =~ \.vcf\.gz$ ]]; then
                            /bin/zcat "$INPUT_VCF" | /bin/sed -E 's/CIPOS=\.,[0-9]+;CIEND=\.,[0-9]+/CIPOS=0,0;CIEND=0,0/g' > fixed_input.vcf
                        else
                            /bin/sed -E 's/CIPOS=\.,[0-9]+;CIEND=\.,[0-9]+/CIPOS=0,0;CIEND=0,0/g' "$INPUT_VCF" > fixed_input.vcf
                        fi
                        #svtools vcftobedpe can take either .vcf or .vcf.gz
                        /opt/conda/envs/python2/bin/svtools vcftobedpe -i fixed_input.vcf -o tmp.bedpe
                        /bin/grep '^#' tmp.bedpe > tmp.header
                        /usr/local/bin/bedtools pairtopair -is -slop "$SLOPE" -type notboth -a tmp.bedpe -b "$BL_BEDPE" | /bin/cat tmp.header /dev/stdin | /opt/conda/envs/python2/bin/svtools bedpetovcf -i /dev/stdin | /opt/conda/envs/python2/bin/svtools vcfsort /dev/stdin > "$OUT_BASE.vcf"

                        /opt/htslib/bin/bgzip $OUT_BASE.vcf
                        /usr/bin/tabix -p vcf $OUT_BASE.vcf.gz
                    else # blocklist_bedpe is not passed.
                        /usr/local/bin/bedtools sort -header -i "$INPUT_VCF" > $OUT_BASE.vcf
                        /opt/htslib/bin/bgzip $OUT_BASE.vcf
                        /usr/bin/tabix -p vcf $OUT_BASE.vcf.gz
                    fi
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: mgibio/basespace_chromoseq:v12
          - class: ResourceRequirement
            ramMin: 7629.39453125
            outdirMin: 10240
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: survivorAnnotateVariants
    in:
      - id: genome_build
        source: genome_build
      - id: input_vcf
        source: filterBlocklistSurvivor/filtered_sv_vcf
      - id: output_tsv_name
        default: SURVIVOR-merged-AnnotSV/tsv
      - id: snps_vcf
        source: snps_vcf
        pickValue: all_non_null
        valueFrom: $([self])
    out:
      - id: sv_variants_tsv
    run:
        class: CommandLineTool
        id: annotsv
        inputs:
          - id: genome_build
            type: string
          - id: input_vcf
            type: File
          - id: output_tsv_name
            default: AnnotSV.tsv
            type: string
          - id: snps_vcf
            type:
                items: File
                type: array
        outputs:
          - id: sv_variants_tsv
            type: File
            outputBinding:
                glob: $(inputs.output_tsv_name)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    /opt/AnnotSV_2.1/bin/AnnotSV -bedtools /usr/bin/bedtools -outputDir "$PWD" \
                    -genomeBuild $(inputs.genome_build) \
                    -SVinputFile $(inputs.input_vcf.path) \
                    -outputFile $(inputs.output_tsv_name) \
                    -vcfFiles $(inputs.snps_vcf.map(function(el) {return el.path}).join(","))
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: mgibio/annotsv-cwl:2.1
          - class: ResourceRequirement
            ramMin: 7629.39453125
            outdirMin: $((10 + Math.round((function(size_of=0){inputs.snps_vcf.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3 + (function(size_of=0){inputs.input_vcf.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3)) * 1024)
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: bcftoolsMergeSvVcfs
    in:
      - id: force_merge
        source: bcftoolsMergeSvVcfs.force_merge
      - id: merge_method
        default: none
      - id: missing_ref
        source: bcftoolsMergeSvVcfs.missing_ref
      - id: output_type
        default: v
      - id: output_vcf_name
        default: bcftools-sv-merged/vcf
      - id: vcfs
        source: sv_vcfs
    out:
      - id: merged_sv_vcf
    run:
        class: CommandLineTool
        id: bcftoolsMerge
        inputs:
          - id: force_merge
            default: true
            type: boolean
          - id: merge_method
            default: none
            type: string
          - id: missing_ref
            default: false
            type: boolean
          - id: output_type
            default: z
            type: string
          - id: output_vcf_name
            default: bcftools_merged.vcf.gz
            type: string
          - id: vcfs
            type:
                items: File
                type: array
        outputs:
          - id: merged_sv_vcf
            type: File
            outputBinding:
                glob: $(inputs.output_vcf_name)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    /opt/bcftools/bin/bcftools merge \
                    --force-samples $(inputs.force_merge) \
                    --merge $(inputs.merge_method) \
                    --missing-to-ref $(inputs.missing_ref) \
                    --output-type $(inputs.output_type) \
                    --output $(inputs.output_vcf_name) \
                    $(inputs.vcfs.map(function(el) {return el.path}).join(" "))
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: mgibio/bcftools-cwl:1.12
          - class: ResourceRequirement
            ramMin: 3814.697265625
            outdirMin: $((10 + Math.round(2 * (function(size_of=0){inputs.vcfs.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3)) * 1024)
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: filterBlocklistBcftools
    in:
      - id: input_vcf
        source: bcftoolsMergeSvVcfs/merged_sv_vcf
      - id: blocklist_bedpe
        source: blocklist_bedpe
      - id: slope
        source: filterBlocklistBcftools.slope
      - id: output_vcf_basename
        default: bcftools-sv-merged
    out:
      - id: filtered_sv_vcf
      - id: filtered_sv_vcf_tbi
    run:
        class: CommandLineTool
        id: filterSvVcfBlocklistBedpe
        inputs:
          - id: input_vcf
            type: File
          - id: blocklist_bedpe
            type:
              - File
              - 'null'
          - id: slope
            default: 100
            type: int
          - id: output_vcf_basename
            default: blocklist_filtered
            type: string
        outputs:
          - id: filtered_sv_vcf
            type: File
            outputBinding:
                glob: $(inputs.output_vcf_basename + ".vcf.gz")
          - id: filtered_sv_vcf_tbi
            type: File
            outputBinding:
                glob: $(inputs.output_vcf_basename + ".vcf.gz.tbi")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eou pipefail
                    set -o errexit

                    INPUT_VCF="$(inputs.input_vcf.path)"
                    SLOPE="$(inputs.slope)"
                    OUT_BASE="$(inputs.output_vcf_basename)"
                    if [[ $(inputs.blocklist_bedpe) = true ]];then # blocklist_bedpe is passed.
                        BL_BEDPE="$(inputs.blocklist_bedpe === null ? "" : inputs.blocklist_bedpe.path)"

                        #CNVkit outputs invalid format like CIPOS=.,894;CIEND=.,894, which can cause svtools vcftobedpe fail
                        if [[ "$INPUT_VCF" =~ \.vcf\.gz$ ]]; then
                            /bin/zcat "$INPUT_VCF" | /bin/sed -E 's/CIPOS=\.,[0-9]+;CIEND=\.,[0-9]+/CIPOS=0,0;CIEND=0,0/g' > fixed_input.vcf
                        else
                            /bin/sed -E 's/CIPOS=\.,[0-9]+;CIEND=\.,[0-9]+/CIPOS=0,0;CIEND=0,0/g' "$INPUT_VCF" > fixed_input.vcf
                        fi
                        #svtools vcftobedpe can take either .vcf or .vcf.gz
                        /opt/conda/envs/python2/bin/svtools vcftobedpe -i fixed_input.vcf -o tmp.bedpe
                        /bin/grep '^#' tmp.bedpe > tmp.header
                        /usr/local/bin/bedtools pairtopair -is -slop "$SLOPE" -type notboth -a tmp.bedpe -b "$BL_BEDPE" | /bin/cat tmp.header /dev/stdin | /opt/conda/envs/python2/bin/svtools bedpetovcf -i /dev/stdin | /opt/conda/envs/python2/bin/svtools vcfsort /dev/stdin > "$OUT_BASE.vcf"

                        /opt/htslib/bin/bgzip $OUT_BASE.vcf
                        /usr/bin/tabix -p vcf $OUT_BASE.vcf.gz
                    else # blocklist_bedpe is not passed.
                        /usr/local/bin/bedtools sort -header -i "$INPUT_VCF" > $OUT_BASE.vcf
                        /opt/htslib/bin/bgzip $OUT_BASE.vcf
                        /usr/bin/tabix -p vcf $OUT_BASE.vcf.gz
                    fi
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: mgibio/basespace_chromoseq:v12
          - class: ResourceRequirement
            ramMin: 7629.39453125
            outdirMin: 10240
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: bcftoolsAnnotateVariants
    in:
      - id: genome_build
        source: genome_build
      - id: input_vcf
        source: filterBlocklistBcftools/filtered_sv_vcf
      - id: output_tsv_name
        default: bcftools-merged-AnnotSV/tsv
      - id: snps_vcf
        source: snps_vcf
        pickValue: all_non_null
        valueFrom: $([self])
    out:
      - id: sv_variants_tsv
    run:
        class: CommandLineTool
        id: annotsv
        inputs:
          - id: genome_build
            type: string
          - id: input_vcf
            type: File
          - id: output_tsv_name
            default: AnnotSV.tsv
            type: string
          - id: snps_vcf
            type:
                items: File
                type: array
        outputs:
          - id: sv_variants_tsv
            type: File
            outputBinding:
                glob: $(inputs.output_tsv_name)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    /opt/AnnotSV_2.1/bin/AnnotSV -bedtools /usr/bin/bedtools -outputDir "$PWD" \
                    -genomeBuild $(inputs.genome_build) \
                    -SVinputFile $(inputs.input_vcf.path) \
                    -outputFile $(inputs.output_tsv_name) \
                    -vcfFiles $(inputs.snps_vcf.map(function(el) {return el.path}).join(","))
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: mgibio/annotsv-cwl:2.1
          - class: ResourceRequirement
            ramMin: 7629.39453125
            outdirMin: $((10 + Math.round((function(size_of=0){inputs.snps_vcf.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3 + (function(size_of=0){inputs.input_vcf.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3)) * 1024)
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: bcftoolsAnnotsvFilter
    in:
      - id: all_CDS
        source: bcftoolsAnnotsvFilter.all_CDS
      - id: annotsv_tsv
        source: bcftoolsAnnotateVariants/sv_variants_tsv
      - id: filtering_frequency
        default: 0/05
      - id: ignore_pass_filter
        source: bcftoolsAnnotsvFilter.ignore_pass_filter
      - id: output_tsv_name
        source: bcftoolsAnnotsvFilter.output_tsv_name
    out:
      - id: filtered_tsv
    run:
        class: CommandLineTool
        id: annotsvFilter
        inputs:
          - id: all_CDS
            default: false
            type: boolean
          - id: annotsv_tsv
            type: File
          - id: filtering_frequency
            default: 0.05
            type: float
          - id: ignore_pass_filter
            default: false
            type: boolean
          - id: output_tsv_name
            default: filtered-bcftools-merged-AnnotSV.tsv
            type: string
        outputs:
          - id: filtered_tsv
            type: File
            outputBinding:
                glob: $(inputs.output_tsv_name)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    python -c '
                    import csv
                    import sys
                    input_file_name  = "$(inputs.annotsv_tsv.path)"
                    output_file_name = "$(inputs.output_tsv_name)"
                    filtering_frequency = $(inputs.filtering_frequency)
                    all_cds = $(inputs.all_CDS ? "True" : "False")
                    ignore_pass_filter = $(inputs.ignore_pass_filter ? "True" : "False")
                    with open(input_file_name, "r") as file_in, open(output_file_name, "w") as file_out:
                        file_in = csv.DictReader(file_in, delimiter="\t")
                        file_out = csv.DictWriter(file_out, fieldnames=file_in.fieldnames, delimiter="\t")
                        file_out.writeheader()
                        total_sv_count = 0
                        pass_sv_count = 0
                        for row in file_in:
                            total_sv_count += 1
                            if(row["AnnotSV type"] == "split" \
                                and (row["FILTER"] == "PASS" or ignore_pass_filter) \
                                and (int(row["CDS length"]) > 0 or all_cds) \
                                and float(row["IMH_AF"]) < filtering_frequency
                                and float(row["1000g_max_AF"]) < filtering_frequency
                                and not(float(row["DGV_LOSS_Frequency"]) > filtering_frequency and "DEL" in row["SV type"])
                                and not(float(row["DGV_GAIN_Frequency"]) < filtering_frequency and ("DUP" in row["SV type"] or "INS" in row["SV type"]))
                                and not(("Manta" in row["ID"] and "IMPRECISE" in row["INFO"]) or (row["QUAL"] != "." and "IMPRECISE" in row["INFO"])) ):
                                file_out.writerow(row)
                                pass_sv_count += 1
                        print("total sv count:",total_sv_count)
                        print("total sv passed count:",pass_sv_count)
                    '
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: python:3
          - class: ResourceRequirement
            ramMin: 3814.697265625
            outdirMin: $((10 + Math.round(2 * (function(size_of=0){inputs.annotsv_tsv.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3)) * 1024)
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
