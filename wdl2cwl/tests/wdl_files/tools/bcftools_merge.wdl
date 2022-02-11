version 1.0

task bcftoolsMerge {
  input {
    Boolean force_merge = true
    String merge_method = "none"  # enum ["none", "snps", "indels", "both", "all", "id"]
    Boolean missing_ref = false
    String output_type = "z"  # enum ["b", "u", "z", "v"]
    String output_vcf_name = "bcftools_merged.vcf.gz"
    Array[File] vcfs
  }

  Int space_needed_gb = 10 + round(2 * size(vcfs, "GB"))
  runtime {
    memory: "4GB"
    docker: "mgibio/bcftools-cwl:1.12"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /opt/bcftools/bin/bcftools merge \
    --force-samples ~{force_merge} \
    --merge ~{merge_method} \
    --missing-to-ref ~{missing_ref} \
    --output-type ~{output_type} \
    --output ~{output_vcf_name} \
    ~{sep=" " vcfs}
  >>>

  output {
    File merged_sv_vcf = output_vcf_name
  }
}
