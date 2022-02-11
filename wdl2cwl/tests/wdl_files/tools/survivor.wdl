version 1.0

task survivor {
  input {
    Array[File] vcfs
    Int max_distance_to_merge
    Int minimum_sv_calls
    Boolean same_type
    Boolean same_strand
    Boolean estimate_sv_distance
    Int minimum_sv_size
    String cohort_name = "SURVIVOR-sv-merged.vcf"
  }

  Int space_needed_gb = 10
  runtime {
    docker: "mgibio/survivor-cwl:1.0.6.2"
    memory: "2GB"
    cpu: 1
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /bin/bash /usr/bin/survivor_merge_helper.sh \
    ~{sep="," vcfs} ~{max_distance_to_merge} ~{minimum_sv_calls} \
    ~{true="1" false="0" same_type} \
    ~{true="1" false="0" same_strand} \
    ~{true="1" false="0" estimate_sv_distance} \
    ~{minimum_sv_size} ~{cohort_name}
  >>>

  output {
    File merged_vcf = cohort_name
  }
}
