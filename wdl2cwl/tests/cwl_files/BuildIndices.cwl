cwlVersion: v1.2
id: BuildIndices
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
inputs:
  - id: gtf_version
    doc: the actual number of gencode, ex.  27
    type: string
  - id: organism
    doc: Either 'human' or 'mouse'
    type: string
  - id: organism_prefix
    doc: Either 'h' or 'm'
    type: string
  - id: genome_short_string
    doc: e.g. hg38, mm10
    type: string
  - id: dbsnp_version
    doc: integer num, ex 150
    type: string
steps:
  - id: _GetReferences.references.genome_fa
    in:
      - id: target
        source: GetReferences/references
        valueFrom: self.genome_fa
    out:
      - genome_fa
    run:
        class: ExpressionTool
        inputs:
          - id: target
            type: Any
        outputs:
          - id: genome_fa
            type: File
        expression: '${return {"genome_fa": self}; }'
  - id: _GetReferences.references.annotation_gtf
    in:
      - id: target
        source: GetReferences/references
        valueFrom: self.annotation_gtf
    out:
      - annotation_gtf
    run:
        class: ExpressionTool
        inputs:
          - id: target
            type: Any
        outputs:
          - id: annotation_gtf
            type: File
        expression: '${return {"annotation_gtf": self}; }'
  - id: _BuildStarSingleNucleus.modified_references.genome_fa
    in:
      - id: target
        source: BuildStarSingleNucleus/modified_references
        valueFrom: self.genome_fa
    out:
      - genome_fa
    run:
        class: ExpressionTool
        inputs:
          - id: target
            type: Any
        outputs:
          - id: genome_fa
            type: File
        expression: '${return {"genome_fa": self}; }'
  - id: _BuildStarSingleNucleus.modified_references.annotation_gtf
    in:
      - id: target
        source: BuildStarSingleNucleus/modified_references
        valueFrom: self.annotation_gtf
    out:
      - annotation_gtf
    run:
        class: ExpressionTool
        inputs:
          - id: target
            type: Any
        outputs:
          - id: annotation_gtf
            type: File
        expression: '${return {"annotation_gtf": self}; }'
  - id: GetReferences
    in:
      - id: gtf_version
        source: gtf_version
      - id: organism
        source: organism
      - id: organism_prefix
        source: organism_prefix
    out:
      - id: references
    run:
        class: CommandLineTool
        id: GetReferences
        doc: Download files needed for building the designated references
        inputs:
          - id: gtf_version
            type: string
          - id: organism
            type: string
          - id: organism_prefix
            type: string
        outputs:
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
            outputBinding:
                outputEval: '$({ "genome_fa": { "class": "File", "path": runtime.outdir+"/"+"GRC"
                    + inputs.organism_prefix + "38.primary_assembly.genome.fa" },
                    "annotation_gtf": { "class": "File", "path": runtime.outdir+"/"+"gencode.v"
                    + inputs.gtf_version + ".primary_assembly.annotation.gtf" } })'
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eo pipefail

                    ## download fasta
                    wget $("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_" + inputs.organism + "/release_" + inputs.gtf_version)/$("GRC" + inputs.organism_prefix + "38.primary_assembly.genome.fa").gz
                    gunzip $("GRC" + inputs.organism_prefix + "38.primary_assembly.genome.fa").gz

                    ## download gtf file
                    wget $("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_" + inputs.organism + "/release_" + inputs.gtf_version)/$("gencode.v" + inputs.gtf_version + ".primary_assembly.annotation.gtf").gz
                    gunzip $("gencode.v" + inputs.gtf_version + ".primary_assembly.annotation.gtf").gz
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/secondary-analysis-star:v2.7.9a
          - class: ResourceRequirement
            outdirMin: 10240
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildPicardRefFlat
    in:
      - id: references
        source: GetReferences/references
    out:
      - id: refflat
    run:
        class: CommandLineTool
        id: BuildPicardRefFlat
        inputs:
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
        outputs:
          - id: refflat
            type: File
            outputBinding:
                glob: $(inputs.references.annotation_gtf.basename.replace(/\.gtf$/,
                    '')  + ".refflat.txt")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4+

                    set -eo pipefail

                    gtfToGenePred -genePredExt -geneNameAsName2  $(inputs.references.annotation_gtf.path) refflat.tmp.txt

                    paste <(cut -f 12 refflat.tmp.txt) <(cut -f 1-10 refflat.tmp.txt) > $(inputs.references.annotation_gtf.basename.replace(/\.gtf$/, '')  + ".refflat.txt")

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/gtf_to_genepred:v0.0.0
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 8192.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildIntervalList
    in:
      - id: references
        source: GetReferences/references
    out:
      - id: interval_list
    run:
        class: CommandLineTool
        id: BuildIntervalList
        inputs:
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
        outputs:
          - id: interval_list
            type: File
            outputBinding:
                glob: $(inputs.references.annotation_gtf.basename.replace(/\.gtf$/,
                    '')  + ".interval_list")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4+

                    set -eo pipefail


                    # index the fasta file

                    mkdir genome_fa
                    ln -s $(inputs.references.genome_fa.path) genome_fa/$(inputs.references.genome_fa.basename)
                    samtools faidx genome_fa/$(inputs.references.genome_fa.basename)
                    cut -f1,2 genome_fa/$(inputs.references.genome_fa.basename).fai > sizes.genome

                    awk -F '\t'  '{  printf "@SQ\tSN:%s\tLN:%s\n", $1, $2 }' sizes.genome  >> $(inputs.references.annotation_gtf.basename.replace(/\.gtf$/, '')  + ".interval_list")

                    grep 'gene_type "rRNA"' $(inputs.references.annotation_gtf.path) |
                        awk '$3 == "transcript"' |
                    cut -f1,4,5,7,9 |
                    perl -lane '
                        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
                        print join "\t", (@F[0,1,2,3], $1)
                    ' |
                    sort -k1V -k2n -k3n  >> $(inputs.references.annotation_gtf.basename.replace(/\.gtf$/, '')  + ".interval_list")

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/secondary-analysis-umitools:0.0.1
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 8192.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildStar
    in:
      - id: gtf_version
        source: gtf_version
      - id: organism
        source: organism
      - id: references
        source: GetReferences/references
    out:
      - id: star_index
    run:
        class: CommandLineTool
        id: BuildStar
        doc: build reference index files for STAR aligner
        inputs:
          - id: gtf_version
            type: string
          - id: organism
            type: string
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
        outputs:
          - id: star_index
            type: File
            outputBinding:
                glob: $("star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version
                    + ".tar")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eo pipefail

                    mkdir star
                    STAR --runMode genomeGenerate \
                      --genomeDir star \
                      --genomeFastaFiles $(inputs.references.genome_fa.path) \
                      --sjdbGTFfile $(inputs.references.annotation_gtf.path) \
                      --sjdbOverhang 100 \
                      --runThreadN 16

                    tar -cvf $("star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version + ".tar") star
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/secondary-analysis-star:v2.7.9a
          - class: ResourceRequirement
            coresMin: 16
            ramMin: 51200.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildStarSingleNucleus
    in:
      - id: gtf_version
        source: gtf_version
      - id: organism
        source: organism
      - id: organism_prefix
        source: organism_prefix
      - id: references
        source: GetReferences/references
    out:
      - id: star_index
      - id: annotation_gtf_modified_introns
      - id: modified_references
    run:
        class: CommandLineTool
        id: BuildStarSingleNucleus
        doc: Modify gtf files and build reference index files for STAR aligner
        inputs:
          - id: gtf_version
            type: string
          - id: organism
            type: string
          - id: organism_prefix
            type: string
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
        outputs:
          - id: star_index
            type: File
            outputBinding:
                glob: $("modified_" + "star_primary_gencode_" + inputs.organism +
                    "_v" + inputs.gtf_version + ".tar")
          - id: annotation_gtf_modified_introns
            type: File
            outputBinding:
                glob: $("introns_modified_gencode.v" + inputs.gtf_version + ".primary_assembly.annotation.gtf")
          - id: modified_references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
            outputBinding:
                outputEval: '$({ "genome_fa": { "class": "File", "path": runtime.outdir+"/"+"modified_GRC"
                    + inputs.organism_prefix + "38.primary_assembly.genome.fa" },
                    "annotation_gtf": { "class": "File", "path": runtime.outdir+"/"+"modified_gencode.v"
                    + inputs.gtf_version + ".primary_assembly.annotation.gtf" } })'
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eo pipefail

                    /script/modify_gtf_$(inputs.organism).sh $(inputs.references.genome_fa.path) $(inputs.references.annotation_gtf.path)

                    mkdir star
                    STAR --runMode genomeGenerate \
                    --genomeDir star \
                    --genomeFastaFiles $("modified_GRC" + inputs.organism_prefix + "38.primary_assembly.genome.fa") \
                    --sjdbGTFfile $("modified_gencode.v" + inputs.gtf_version + ".primary_assembly.annotation.gtf") \
                    --sjdbOverhang 100 \
                    --runThreadN 16

                    tar -cvf $("modified_" + "star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version + ".tar") star

                    python3  /script/add-introns-to-gtf.py   --input-gtf $("modified_gencode.v" + inputs.gtf_version + ".primary_assembly.annotation.gtf")  --output-gtf $("introns_modified_gencode.v" + inputs.gtf_version + ".primary_assembly.annotation.gtf")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: 'quay.io/humancellatlas/snss2-indices:1.1.0 '
          - class: ResourceRequirement
            coresMin: 16
            ramMin: 51200.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildRsem
    in:
      - id: gtf_version
        source: gtf_version
      - id: organism
        source: organism
      - id: references
        source: GetReferences/references
    out:
      - id: rsem_index
    run:
        class: CommandLineTool
        id: BuildRsem
        doc: build reference index files for RSEM
        inputs:
          - id: gtf_version
            type: string
          - id: organism
            type: string
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
        outputs:
          - id: rsem_index
            type: File
            outputBinding:
                glob: $("rsem_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version
                    + ".tar")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eo pipefail
                    mkdir rsem
                    rsem-prepare-reference --gtf $(inputs.references.annotation_gtf.path) --bowtie $(inputs.references.genome_fa.path) rsem/rsem_trans_index
                    tar -cvf $("rsem_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version + ".tar") rsem/
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0
          - class: ResourceRequirement
            ramMin: 10240.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildHisat2FromRsem
    in:
      - id: rsem_index
        source: BuildRsem/rsem_index
    out:
      - id: hisat2_index
    run:
        class: CommandLineTool
        id: BuildHisat2FromRsem
        inputs:
          - id: rsem_index
            type: File
        outputs:
          - id: hisat2_index
            type: File
            outputBinding:
                glob: $("hisat2_from_rsem_" + inputs.rsem_index.basename.replace(/\.tar$/,
                    '')  + ".tar.gz")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4


                    # extract rsem index
                    tar -xf $(inputs.rsem_index.path)

                    # build index
                    hisat2-build -p 8 rsem/rsem_trans_index.idx.fa $("hisat2_from_rsem_" + inputs.rsem_index.basename.replace(/\.tar$/, '') )
                    mkdir $("hisat2_from_rsem_" + inputs.rsem_index.basename.replace(/\.tar$/, '') )
                    mv ./*.ht2 $("hisat2_from_rsem_" + inputs.rsem_index.basename.replace(/\.tar$/, '') )
                    tar -zcvf $("hisat2_from_rsem_" + inputs.rsem_index.basename.replace(/\.tar$/, '')  + ".tar.gz") $("hisat2_from_rsem_" + inputs.rsem_index.basename.replace(/\.tar$/, '') )
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 8192.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildHisat2
    in:
      - id: gtf_version
        source: gtf_version
      - id: organism
        source: organism
      - id: genome_fa
        source: GetReferences/references
        valueFrom: $(self.genome_fa)
    out:
      - id: hisat2_index
    run:
        class: CommandLineTool
        id: BuildHisat2
        inputs:
          - id: gtf_version
            type: string
          - id: organism
            type: string
          - id: genome_fa
            type: File
        outputs:
          - id: hisat2_index
            type: File
            outputBinding:
                glob: $("hisat2_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version
                    + ".tar.gz")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eo pipefail

                    # build index
                    hisat2-build -p 8 $(inputs.genome_fa.path) $("hisat2_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version)
                    mkdir $("hisat2_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version)
                    mv ./*.ht2 $("hisat2_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version)
                    tar -zcvf $("hisat2_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version + ".tar.gz") $("hisat2_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version)
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 65536.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: BuildHisat2SnpHaplotypeSplicing
    in:
      - id: organism
        source: organism
      - id: genome_short_string
        source: genome_short_string
      - id: references
        source: GetReferences/references
      - id: gtf_version
        source: gtf_version
      - id: dbsnp_version
        source: dbsnp_version
    out:
      - id: hisat2_index
    run:
        class: CommandLineTool
        id: BuildHisat2SnpHaplotypeSplicing
        inputs:
          - id: organism
            type: string
          - id: genome_short_string
            type: string
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
          - id: gtf_version
            type: string
          - id: dbsnp_version
            type: string
        outputs:
          - id: hisat2_index
            type: File
            outputBinding:
                glob: $("star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version
                    + ".tar.gz")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4+


                    HISAT2_DIR=/opt/tools/hisat2-2.1.0

                    # Compressed fasta required here
                    gzip $(inputs.references.genome_fa.path)

                    # download snp file
                    wget http://hgdownload.cse.ucsc.edu/goldenPath/$(inputs.genome_short_string)/database/$("snp" + inputs.dbsnp_version + "Common.txt").gz
                    gunzip $("snp" + inputs.dbsnp_version + "Common.txt").gz

                    # extract snps, splice sites, and exon information
                    $HISAT2_DIR/hisat2_extract_snps_UCSC.py $(inputs.references.genome_fa.path).gz $("snp" + inputs.dbsnp_version + "Common.txt") genome
                    $HISAT2_DIR/hisat2_extract_splice_sites.py $(inputs.references.annotation_gtf.path) > genome.ss
                    $HISAT2_DIR/hisat2_extract_exons.py $(inputs.references.annotation_gtf.path) > genome.exon

                    # build the hisat2 reference
                    $HISAT2_DIR/hisat2-build \
                      -p 8 \
                      genome.fa \
                      --snp genome.snp \
                      --haplotype genome.haplotype \
                      --ss genome.ss \
                      --exon genome.exon \
                      genome_snp_tran

                    mkdir $("star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version)
                    cp ./*.ht2 $("star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version)
                    tar -zcvf $("star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version + ".tar.gz") $("star_primary_gencode_" + inputs.organism + "_v" + inputs.gtf_version)

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/humancellatlas/secondary-analysis-hisat2:v0.3.0-2-2.1.0
          - class: ResourceRequirement
            coresMin: 16
            ramMin: 245760.0
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
outputs:
  - id: BuildIndices.star_index
    outputSource: BuildStar/star_index
    type: File
  - id: BuildIndices.snSS2_star_index
    outputSource: BuildStarSingleNucleus/star_index
    type: File
  - id: BuildIndices.rsem_index
    outputSource: BuildRsem/rsem_index
    type: File
  - id: BuildIndices.hisat2_from_rsem_index
    outputSource: BuildHisat2FromRsem/hisat2_index
    type: File
  - id: BuildIndices.hisat2_index
    outputSource: BuildHisat2/hisat2_index
    type: File
  - id: BuildIndices.hisat2_snp_haplotype_splicing_index
    outputSource: BuildHisat2SnpHaplotypeSplicing/hisat2_index
    type: File
  - id: BuildIndices.refflat
    outputSource: BuildPicardRefFlat/refflat
    type: File
  - id: BuildIndices.interval_list
    outputSource: BuildIntervalList/interval_list
    type: File
  - id: BuildIndices.genome_fa
    outputSource: _GetReferences.references.genome_fa/genome_fa
    type: File
  - id: BuildIndices.annotation_gtf
    outputSource: _GetReferences.references.annotation_gtf/annotation_gtf
    type: File
  - id: BuildIndices.snSS2_genome_fa
    outputSource: _BuildStarSingleNucleus.modified_references.genome_fa/genome_fa
    type: File
  - id: BuildIndices.snSS2_annotation_gtf
    outputSource: _BuildStarSingleNucleus.modified_references.annotation_gtf/annotation_gtf
    type: File
  - id: BuildIndices.snSS2_annotation_gtf_introns
    outputSource: BuildStarSingleNucleus/annotation_gtf_modified_introns
    type: File
