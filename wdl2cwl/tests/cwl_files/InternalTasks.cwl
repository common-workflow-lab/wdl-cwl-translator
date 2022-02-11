cwlVersion: v1.2
$graph:
  - class: CommandLineTool
    id: MakeSafeFilename
    inputs:
      - id: name
        type: string
    outputs:
      - id: output_safe_name
        type: string
        outputBinding:
            loadContents: true
            glob: safe_name.txt
            outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                echo '$(inputs.name)' | tr ' "#$%&*/:;<=>?@[]^{}|~\\()' '_' > safe_name.txt
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: gcr.io/gcp-runtimes/ubuntu_16_0_4:latest
      - class: ResourceRequirement
        ramMin: 1024.0
        outdirMin: 10240
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: DownloadGenotypes
    inputs:
      - id: sample_alias
        type: string
      - id: sample_lsid
        type: string
      - id: output_vcf_base_name
        type: string
      - id: compress
        default: true
        type: boolean
      - id: ignoreSpecificGenotypesLsid
        type:
          - string
          - 'null'
      - id: ignoreSpecificGenotypesPlatform
        type:
          - string
          - 'null'
      - id: haplotype_database_file
        type: File
      - id: ref_fasta
        type: File
      - id: ref_fasta_index
        type: File
      - id: ref_dict
        type: File
      - id: environment
        type: string
      - id: vault_token_path
        type: File
      - id: max_retries
        type:
          - int
          - 'null'
      - id: preemptible_tries
        type:
          - int
          - 'null'
    outputs:
      - id: fingerprint_retrieved
        type: boolean
        outputBinding:
            loadContents: true
            glob: $("fp_retrieved.txt")
            outputEval: |-
                ${
                var contents = self[0].contents.trim().toLowerCase()
                if (contents == 'true') { return true;}
                if (contents == 'false') { return false;}
                throw "'read_boolean' received neither 'true' nor 'false': " + self[0].contents;
                }
      - id: reference_fingerprint_vcf
        type: File
        outputBinding:
            glob: '$(inputs.output_vcf_base_name + inputs.compress ? ".vcf.gz" : ".vcf")'
      - id: reference_fingerprint_vcf_index
        type: File
        outputBinding:
            glob: '$(inputs.output_vcf_base_name + inputs.compress ? ".vcf.gz" : ".vcf"
                + inputs.compress ? ".tbi" : ".idx")'
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4


                export VAULT_ADDR=https://clotho.broadinstitute.org:8200
                export VAULT_TOKEN=\$(cat $(inputs.vault_token_path.path))
                if [ $(inputs.environment) == prod ]; then
                  export MERCURY_AUTH_KEY=secret/dsde/gotc/prod/wdl/secrets
                  export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal/mercury-ws/fingerprint
                else
                  export MERCURY_AUTH_KEY=secret/dsde/gotc/dev/wdl/secrets
                  export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal-test/mercury-ws/fingerprint
                fi

                exit_code=0

                # TODO - there is a bug in DownloadGenotypes - we should NOT have to set EXPECTED_GENOTYPING_PLATFORMS here
                # it * should * default to all of them.

                java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
                DownloadGenotypes \
                  --SAMPLE_ALIAS "$(inputs.sample_alias)" \
                  --SAMPLE_LSID "$(inputs.sample_lsid)" \
                  --OUTPUT "$(inputs.output_vcf_base_name + inputs.compress ? ".vcf.gz" : ".vcf")" \
                  --CREATE_INDEX true \
                  --REFERENCE_SEQUENCE $(inputs.ref_fasta.path) \
                  --HAPLOTYPE_MAP $(inputs.haplotype_database_file.path) \
                  --EXPECTED_GENOTYPING_PLATFORMS FLUIDIGM \
                  --EXPECTED_GENOTYPING_PLATFORMS GENERAL_ARRAY \
                  $(inputs.ignoreSpecificGenotypesLsid ? inputs.ignoreSpecificGenotypesLsid === null ? "" : "--IGNORE_SPECIFIC_GENOTYPES_LSID "" + inputs.ignoreSpecificGenotypesLsid + """ : "") \
                  $(inputs.ignoreSpecificGenotypesPlatform ? inputs.ignoreSpecificGenotypesPlatform === null ? "" : "--IGNORE_SPECIFIC_GENOTYPES_PLATFORM "" + inputs.ignoreSpecificGenotypesPlatform + """ : "") \
                  --MERCURY_FP_STORE_URI $MERCURY_FP_STORE_URI \
                  --CREDENTIALS_VAULT_PATH $MERCURY_AUTH_KEY \
                  --ERR_NO_GENOTYPES_AVAILABLE 7
                exit_code=$?

                if [ $exit_code -eq 0 ]; then
                  echo "true" > $("fp_retrieved.txt")
                elif [ $exit_code -eq 7 ]; then
                  # Exit code from DownloadGenotypes if no fingerprints were found.
                  # Treat this as a normal condition, but set a variable to indicate no fingerprints available.
                  # Create empty file so that it exists.
                  exit_code=0
                  echo "Found no fingerprints for $(inputs.sample_lsid)"
                  echo "false" > $("fp_retrieved.txt")
                  touch $(inputs.output_vcf_base_name + inputs.compress ? ".vcf.gz" : ".vcf")
                  touch $(inputs.output_vcf_base_name + inputs.compress ? ".vcf.gz" : ".vcf" + inputs.compress ? ".tbi" : ".idx")
                else
                  echo "false" > $("fp_retrieved.txt")
                fi

                exit $exit_code
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.0-1641925612
      - class: ResourceRequirement
        ramMin: 3500.0
        outdirMin: 1024
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: UploadFingerprintToMercury
    inputs:
      - id: fingerprint_json_file
        type: File
      - id: gtc_file
        type: File
      - id: environment
        type: string
      - id: vault_token_path
        type: File
      - id: max_retries
        type:
          - int
          - 'null'
      - id: preemptible_tries
        type:
          - int
          - 'null'
    outputs: []
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -eo pipefail

                export VAULT_ADDR=https://clotho.broadinstitute.org:8200
                export VAULT_TOKEN=\$(cat $(inputs.vault_token_path.path))
                if [ $(inputs.environment) == prod ]; then
                  export MERCURY_AUTH_KEY=secret/dsde/gotc/prod/wdl/secrets
                  export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal/mercury-ws/fingerprint
                else
                  export MERCURY_AUTH_KEY=secret/dsde/gotc/dev/wdl/secrets
                  export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal-test/mercury-ws/fingerprint
                fi

                du -k $(inputs.gtc_file.path) | cut -f 1 > size.txt

                # TODO -Fix UploadFingerprintToMercury so I don't need to pass a file size

                java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
                  UploadFingerprintToMercury \
                  --INPUT "$(inputs.fingerprint_json_file.path)" \
                  --GTC_FILE_SIZE size.txt \
                  --MERCURY_FP_STORE_URI $MERCURY_FP_STORE_URI \
                  --CREDENTIALS_VAULT_PATH $MERCURY_AUTH_KEY \
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.0-1641925612
      - class: ResourceRequirement
        ramMin: 3500.0
        outdirMin: 1024
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
