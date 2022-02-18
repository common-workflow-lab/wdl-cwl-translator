cwlVersion: v1.2
id: gwas
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
inputs:
  - id: vcf
    type: File
  - id: metadata_csv
    type: File
steps:
  - id: parse_metadata
    in:
      - id: metadata_csv
        source: metadata_csv
    out:
      - id: covariates
      - id: phenotypes
      - id: sex
      - id: ids
    run:
        class: CommandLineTool
        id: parse_metadata
        inputs:
          - id: metadata_csv
            type: File
        outputs:
          - id: covariates
            type: File
            outputBinding:
                glob: covariates.txt
          - id: phenotypes
            type: File
            outputBinding:
                glob: phenotypes.txt
          - id: sex
            type: File
            outputBinding:
                glob: sex.txt
          - id: ids
            type: File
            outputBinding:
                glob: ids.txt
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    parse_metadata.sh \
                    	-c $(inputs.metadata_csv.path)
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: dnastack/plink:1.9
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3576.2786865234375
            outdirMin: 20480
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: run_gwas
    in:
      - id: vcf
        source: vcf
      - id: covariates
        source: parse_metadata/covariates
      - id: phenotypes
        source: parse_metadata/phenotypes
      - id: sex
        source: parse_metadata/sex
      - id: ids
        source: parse_metadata/ids
    out:
      - id: logistic
    run:
        class: CommandLineTool
        id: run_gwas
        inputs:
          - id: vcf
            type: File
          - id: covariates
            type: File
          - id: phenotypes
            type: File
          - id: sex
            type: File
          - id: ids
            type: File
        outputs:
          - id: logistic
            type: File
            outputBinding:
                glob: $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/,
                    '') + ".assoc.logistic")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    plink \
                    	--vcf $(inputs.vcf.path) \
                    	--maf 0.10 \
                    	--update-ids $(inputs.ids.path) \
                    	--make-bed \
                    	--out $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/, ''))

                    plink \
                    	--bfile $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/, '')) \
                    	--update-sex $(inputs.sex.path) \
                    	--pheno $(inputs.phenotypes.path) \
                    	--make-bed \
                    	--out $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/, ''))

                    # Recode covariates to binary
                    plink \
                    	--bfile $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/, '')) \
                    	--covar $(inputs.covariates.path) \
                    	--dummy-coding \
                    	--write-covar

                    plink \
                    	--bfile $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/, '')) \
                    	--logistic \
                    	--covar plink.cov \
                    	--out $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/, ''))

                    sed -i -e 's/\s\+/,/g' -e 's/^,//g' -e 's/,$//g' $(inputs.vcf.path.replace("\\?.*", "") .split('/').reverse()[0].replace(/\.vcf\.gz$/, '')).assoc.logistic
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: dnastack/plink:1.9
          - class: ResourceRequirement
            coresMin: 4
            ramMin: 15258.7890625
            outdirMin: 358400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: create_plot
    in:
      - id: assoc_file
        source: run_gwas/logistic
    out:
      - id: manhattan_plot
    run:
        class: CommandLineTool
        id: create_plot
        inputs:
          - id: assoc_file
            type: File
        outputs:
          - id: manhattan_plot
            type: File
            outputBinding:
                glob: $(inputs.assoc_file.basename.replace(/\.assoc\.logistic$/, '')  +
                    ".png")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    manhattan_plot.py \
                    	-i $(inputs.assoc_file.path) \
                    	-o $(inputs.assoc_file.basename.replace(/\.assoc\.logistic$/, '') ).png
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: dnastack/plink:1.9
          - class: ResourceRequirement
            coresMin: 2
            ramMin: 7152.557373046875
            outdirMin: 20480
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
outputs:
  - id: gwas.logistic
    outputSource: run_gwas/logistic
    type: File
  - id: gwas.manhattan_plot
    outputSource: create_plot/manhattan_plot
    type: File
