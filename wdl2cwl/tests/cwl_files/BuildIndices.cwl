class: CommandLineTool
id: BuildIntervalList
inputs:
  - id: gtf_version
    type: string
  - id: organism
    type: string
  - id: organism_prefix
    type: string
  - id: gtf_version
    type: string
  - id: organism
    type: string
  - id: references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
  - id: gtf_version
    type: string
  - id: organism
    type: string
  - id: organism_prefix
    type: string
  - id: references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
  - id: gtf_version
    type: string
  - id: organism
    type: string
  - id: references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
  - id: rsem_index
    type: File
  - id: gtf_version
    type: string
  - id: organism
    type: string
  - id: genome_fa
    type: File
  - id: organism
    type: string
  - id: genome_short_string
    type: string
  - id: references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
  - id: gtf_version
    type: string
  - id: dbsnp_version
    type: string
  - id: references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
  - id: references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
outputs:
  - id: references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
  - id: star_index
    type: File
    outputBinding:
        glob: ''
  - id: star_index
    type: File
    outputBinding:
        glob: ''
  - id: annotation_gtf_modified_introns
    type: File
    outputBinding:
        glob: ''
  - id: modified_references
    type:
        fields:
          - name: genome_fa
            type: File
          - name: annotation_gtf
            type: File
        type: References
  - id: rsem_index
    type: File
    outputBinding:
        glob: ''
  - id: hisat2_index
    type: File
    outputBinding:
        glob: ''
  - id: hisat2_index
    type: File
    outputBinding:
        glob: ''
  - id: hisat2_index
    type: File
    outputBinding:
        glob: ''
  - id: refflat
    type: File
    outputBinding:
        glob: ''
  - id: interval_list
    type: File
    outputBinding:
        glob: ''
  - id: pipeline_version
    type: string
    outputBinding:
        glob: 0.1.0
  - id: star_index
    type: File
    outputBinding:
        glob: ''
  - id: snSS2_star_index
    type: File
    outputBinding:
        glob: ''
  - id: rsem_index
    type: File
    outputBinding:
        glob: ''
  - id: hisat2_from_rsem_index
    type: File
    outputBinding:
        glob: ''
  - id: hisat2_index
    type: File
    outputBinding:
        glob: ''
  - id: hisat2_snp_haplotype_splicing_index
    type: File
    outputBinding:
        glob: ''
  - id: refflat
    type: File
    outputBinding:
        glob: ''
  - id: interval_list
    type: File
    outputBinding:
        glob: ''
  - id: genome_fa
    type: File
    outputBinding:
        glob: ''
  - id: annotation_gtf
    type: File
    outputBinding:
        glob: ''
  - id: snSS2_genome_fa
    type: File
    outputBinding:
        glob: ''
  - id: snSS2_annotation_gtf
    type: File
    outputBinding:
        glob: ''
  - id: snSS2_annotation_gtf_introns
    type: File
    outputBinding:
        glob: ''
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/humancellatlas/secondary-analysis-umitools:0.0.1
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4+

            set -eo pipefail


            # index the fasta file

            samtools faidx $(inputs.references.genome_fa)
            cut -f1,2 $(inputs.references.genome_fa).fai > sizes.genome

            awk -F '\t'  '{  printf "@SQ\tSN:%s\tLN:%s\n", $1, $2 }' sizes.genome  >> $(inputs.interval_list_name)

            grep 'gene_type "rRNA"' $(inputs.references.annotation_gtf) |
                awk '$3 == "transcript"' |
            cut -f1,4,5,7,9 |
            perl -lane '
                /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
                print join "\t", (@F[0,1,2,3], $1)
            ' |
            sort -k1V -k2n -k3n  >> $(inputs.interval_list_name)

  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    ramMin: 8192
  - class: ResourceRequirement
    coresMin: ''
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
