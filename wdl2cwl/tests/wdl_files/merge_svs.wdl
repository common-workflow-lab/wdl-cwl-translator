# source https://github.com/griffithlab/analysis-wdls/blob/25cf2010055d79bf7df63ab832363066eba87a26/definitions/subworkflows/merge_svs.wdll

version 1.0

import "tools/survivor.wdl" as s
import "tools/filter_sv_vcf_blocklist_bedpe.wdl" as fsvbb
import "tools/annotsv.wdl" as a
import "tools/bcftools_merge.wdl" as bm
import "tools/annotsv_filter.wdl" as af

workflow mergeSvs {
  input {
    Boolean estimate_sv_distance
    String genome_build
    Int max_distance_to_merge
    Int minimum_sv_calls
    Int minimum_sv_size
    Boolean same_strand
    Boolean same_type
    File? snps_vcf
    Array[File] sv_vcfs
    File? blocklist_bedpe
  }

  call s.survivor as survivorMergeSvVcfs {
    input:
    vcfs=sv_vcfs,
    max_distance_to_merge=max_distance_to_merge,
    minimum_sv_calls=minimum_sv_calls,
    same_type=same_type,
    same_strand=same_strand,
    estimate_sv_distance=estimate_sv_distance,
    minimum_sv_size=minimum_sv_size,
    cohort_name="SURVIVOR-sv-merged.vcf"
  }

  call fsvbb.filterSvVcfBlocklistBedpe as filterBlocklistSurvivor {
    input:
    input_vcf=survivorMergeSvVcfs.merged_vcf,
    blocklist_bedpe=blocklist_bedpe,
    output_vcf_basename="SURVIVOR-sv-merged"
  }

  call a.annotsv as survivorAnnotateVariants {
    input:
    genome_build=genome_build,
    input_vcf=filterBlocklistSurvivor.filtered_sv_vcf,
    output_tsv_name="SURVIVOR-merged-AnnotSV.tsv",
    snps_vcf=select_all([snps_vcf])
  }

  call bm.bcftoolsMerge as bcftoolsMergeSvVcfs {
    input:
    merge_method="none",
    output_type="v",
    output_vcf_name="bcftools-sv-merged.vcf",
    vcfs=sv_vcfs
  }

  call fsvbb.filterSvVcfBlocklistBedpe as filterBlocklistBcftools {
    input:
    input_vcf=bcftoolsMergeSvVcfs.merged_sv_vcf,
    blocklist_bedpe=blocklist_bedpe,
    output_vcf_basename="bcftools-sv-merged"
  }

  call a.annotsv as bcftoolsAnnotateVariants {
    input:
    genome_build=genome_build,
    input_vcf=filterBlocklistBcftools.filtered_sv_vcf,
    output_tsv_name="bcftools-merged-AnnotSV.tsv",
    snps_vcf=select_all([snps_vcf])
  }

  call af.annotsvFilter as bcftoolsAnnotsvFilter {
    input:
    annotsv_tsv=bcftoolsAnnotateVariants.sv_variants_tsv,
    filtering_frequency="0.05"
  }

  output {
    File bcftools_merged_sv_vcf = filterBlocklistBcftools.filtered_sv_vcf
    File bcftools_merged_annotated_tsv = bcftoolsAnnotateVariants.sv_variants_tsv
    File bcftools_merged_filtered_annotated_tsv = bcftoolsAnnotsvFilter.filtered_tsv
    File survivor_merged_sv_vcf = filterBlocklistSurvivor.filtered_sv_vcf
    File survivor_merged_annotated_tsv = survivorAnnotateVariants.sv_variants_tsv
  }
}
