version 1.0

task filterSvVcfBlocklistBedpe {
  input {
    File input_vcf
    File? blocklist_bedpe
    Int slope = 100
    String output_vcf_basename = "blocklist_filtered"
  }

  Int space_needed_gb = 10
  runtime {
    memory: "8GB"
    docker: "mgibio/basespace_chromoseq:v12"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    set -eou pipefail
    set -o errexit

    INPUT_VCF="~{input_vcf}"
    SLOPE="~{slope}"
    OUT_BASE="~{output_vcf_basename}"
    if [[ ~{defined(blocklist_bedpe)} = true ]];then # blocklist_bedpe is passed.
        BL_BEDPE="~{blocklist_bedpe}"

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
  >>>

  output {
    File filtered_sv_vcf = output_vcf_basename + ".vcf.gz"
    File filtered_sv_vcf_tbi = output_vcf_basename + ".vcf.gz.tbi"
  }
}
