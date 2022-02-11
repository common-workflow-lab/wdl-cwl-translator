version 1.0

task annotsv {
  input {
    String genome_build
    File input_vcf
    String output_tsv_name = "AnnotSV.tsv"
    Array[File] snps_vcf
  }

  Int space_needed_gb = 10 + round(size(snps_vcf, "GB") + size(input_vcf, "GB"))
  runtime {
    memory: "8GB"
    docker: "mgibio/annotsv-cwl:2.1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /opt/AnnotSV_2.1/bin/AnnotSV -bedtools /usr/bin/bedtools -outputDir "$PWD" \
    -genomeBuild ~{genome_build} \
    -SVinputFile ~{input_vcf} \
    -outputFile ~{output_tsv_name} \
    -vcfFiles ~{sep="," snps_vcf}
  >>>

  output {
    File sv_variants_tsv = output_tsv_name
  }
}
