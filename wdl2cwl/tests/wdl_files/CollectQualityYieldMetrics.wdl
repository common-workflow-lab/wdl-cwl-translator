version 1.0

## Source: https://github.com/broadinstitute/warp/blob/9ea4dc8ad4eca794ed336451fa8b23367132a820/tasks/broad/Qc.wdl#L19
##
## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for QC of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
  input {
    File input_bam
    String metrics_filename
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20

  command {
    java -Xms2000m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=~{input_bam} \
      OQ=true \
      OUTPUT=~{metrics_filename}
    sed -i -e 1,5d ~{metrics_filename}  # for reproducibility
  }
  
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "40 GiB"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }
  output {
    File quality_yield_metrics = "~{metrics_filename}"
  }
}


