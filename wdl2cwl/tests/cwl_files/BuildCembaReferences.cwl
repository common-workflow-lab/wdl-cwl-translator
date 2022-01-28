class: Workflow
id: BuildCembaReferences
inputs:
  - id: BuildBisulfiteReferences.fasta_input
    type: File
  - id: BuildBisulfiteReferences.monitoring_script
    type:
      - File
      - 'null'
  - id: Bowtie2Build.fasta_input
    type: File
  - id: Bowtie2Build.index_prefix
    type: string
  - id: Bowtie2Build.monitoring_script
    type:
      - File
      - 'null'
  - id: Bowtie2Build.fasta_input
    type: File
  - id: Bowtie2Build.index_prefix
    type: string
  - id: Bowtie2Build.monitoring_script
    type:
      - File
      - 'null'
  - id: CreateReferenceDictionary.reference_fasta
    type: File
  - id: CreateReferenceDictionary.monitoring_script
    type:
      - File
      - 'null'
  - id: CreateReferenceFastaIndex.reference_fasta
    type: File
  - id: CreateReferenceFastaIndex.monitoring_script
    type:
      - File
      - 'null'
outputs:
  - id: fwd_converted_reference_fasta_output
    outputSource: BuildBisulfiteReferences/fwd_converted_reference_fasta_output
    type: File
  - id: rev_converted_reference_fasta_output
    outputSource: BuildBisulfiteReferences/rev_converted_reference_fasta_output
    type: File
  - id: monitoring_log
    outputSource: BuildBisulfiteReferences/monitoring_log
    type: File
  - id: bowtie2_index_files
    outputSource: Bowtie2Build/bowtie2_index_files
    type:
        items: File
        type: array
  - id: monitoring_log
    outputSource: Bowtie2Build/monitoring_log
    type: File
  - id: bowtie2_index_files
    outputSource: Bowtie2Build/bowtie2_index_files
    type:
        items: File
        type: array
  - id: monitoring_log
    outputSource: Bowtie2Build/monitoring_log
    type: File
  - id: ref_dict_output
    outputSource: CreateReferenceDictionary/ref_dict_output
    type: File
  - id: monitoring_log
    outputSource: CreateReferenceDictionary/monitoring_log
    type: File
  - id: ref_index_output
    outputSource: CreateReferenceFastaIndex/ref_index_output
    type: File
  - id: monitoring_log
    outputSource: CreateReferenceFastaIndex/monitoring_log
    type: File
cwlVersion: v1.2
steps:
  - id: BuildBisulfiteReferences
    in:
      - id: fasta_input
        source: BuildBisulfiteReferences.fasta_input
      - id: monitoring_script
        source: BuildBisulfiteReferences.monitoring_script
    out:
      - id: fwd_converted_reference_fasta_output
      - id: rev_converted_reference_fasta_output
      - id: monitoring_log
    run:
        class: CommandLineTool
        id: BuildBisulfiteReferences
        inputs:
          - id: fasta_input
            type: File
          - id: monitoring_script
            type:
              - File
              - 'null'
        outputs:
          - id: fwd_converted_reference_fasta_output
            type: File
            outputBinding:
                glob: $("genome_mfa.CT_conversion.fa")
          - id: rev_converted_reference_fasta_output
            type: File
            outputBinding:
                glob: $("genome_mfa.GA_conversion.fa")
          - id: monitoring_log
            type: File
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/bisulfite-references:1.0
          - class: InitialWorkDirRequirement
            listing:
              - entryname: example.sh
                entry: |4

                    # if the WDL/task contains a monitoring script as input
                    if [ ! -z "$(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)" ]; then
                      chmod a+x $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)
                      $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path) > monitoring.log &
                    else
                      echo "No monitoring script given as input" > monitoring.log &
                    fi

                    python /build_bisulfite_references.py \
                      --input-fasta $(inputs.fasta_input.path) \
                      --forward-convert-out $("genome_mfa.CT_conversion.fa") \
                      --reverse-convert-out $("genome_mfa.GA_conversion.fa")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3337.860107421875
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - example.sh
  - id: Bowtie2Build
    in:
      - id: fasta_input
        source: Bowtie2Build.fasta_input
      - id: index_prefix
        source: Bowtie2Build.index_prefix
      - id: monitoring_script
        source: Bowtie2Build.monitoring_script
    out:
      - id: bowtie2_index_files
      - id: monitoring_log
    run:
        class: CommandLineTool
        id: Bowtie2Build
        inputs:
          - id: fasta_input
            type: File
          - id: index_prefix
            type: string
          - id: monitoring_script
            type:
              - File
              - 'null'
        outputs:
          - id: bowtie2_index_files
            type:
                items: File
                type: array
            outputBinding:
                glob: $(inputs.index_prefix + "*")
          - id: monitoring_log
            type: File
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/bowtie2:2.3.4.3
          - class: InitialWorkDirRequirement
            listing:
              - entryname: example.sh
                entry: |4

                    # if the WDL/task contains a monitoring script as input
                    if [ ! -z "$(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)" ]; then
                      chmod a+x $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)
                      $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path) > monitoring.log &
                    else
                      echo "No monitoring script given as input" > monitoring.log &
                    fi

                    bowtie2-build \
                      -f $(inputs.fasta_input.path) $(inputs.index_prefix)
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 6675.72021484375
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - example.sh
  - id: Bowtie2Build
    in:
      - id: fasta_input
        source: Bowtie2Build.fasta_input
      - id: index_prefix
        source: Bowtie2Build.index_prefix
      - id: monitoring_script
        source: Bowtie2Build.monitoring_script
    out:
      - id: bowtie2_index_files
      - id: monitoring_log
    run:
        class: CommandLineTool
        id: Bowtie2Build
        inputs:
          - id: fasta_input
            type: File
          - id: index_prefix
            type: string
          - id: monitoring_script
            type:
              - File
              - 'null'
        outputs:
          - id: bowtie2_index_files
            type:
                items: File
                type: array
            outputBinding:
                glob: $(inputs.index_prefix + "*")
          - id: monitoring_log
            type: File
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/bowtie2:2.3.4.3
          - class: InitialWorkDirRequirement
            listing:
              - entryname: example.sh
                entry: |4

                    # if the WDL/task contains a monitoring script as input
                    if [ ! -z "$(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)" ]; then
                      chmod a+x $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)
                      $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path) > monitoring.log &
                    else
                      echo "No monitoring script given as input" > monitoring.log &
                    fi

                    bowtie2-build \
                      -f $(inputs.fasta_input.path) $(inputs.index_prefix)
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 6675.72021484375
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - example.sh
  - id: CreateReferenceDictionary
    in:
      - id: reference_fasta
        source: CreateReferenceDictionary.reference_fasta
      - id: monitoring_script
        source: CreateReferenceDictionary.monitoring_script
    out:
      - id: ref_dict_output
      - id: monitoring_log
    run:
        class: CommandLineTool
        id: CreateReferenceDictionary
        inputs:
          - id: reference_fasta
            type: File
          - id: monitoring_script
            type:
              - File
              - 'null'
        outputs:
          - id: ref_dict_output
            type: File
            outputBinding:
                glob: $(inputs.reference_fasta.basename.replace(/\.fa$/, '') .split('/').reverse()[0].replace(/\.fasta$/,
                    '') + ".dict")
          - id: monitoring_log
            type: File
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/picard:2.18.23
          - class: InitialWorkDirRequirement
            listing:
              - entryname: example.sh
                entry: |4+

                    set -euo pipefail

                    # if the WDL/task contains a monitoring script as input
                    if [ ! -z "$(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)" ]; then
                      chmod a+x $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)
                      $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path) > monitoring.log &
                    else
                      echo "No monitoring script given as input" > monitoring.log &
                    fi

                    # create a reference dict
                    java -jar /picard-tools/picard.jar CreateSequenceDictionary \
                      REFERENCE=$(inputs.reference_fasta.path) \
                      OUTPUT=$(inputs.reference_fasta.basename.replace(/\.fa$/, '') .split('/').reverse()[0].replace(/\.fasta$/, '') + ".dict")

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3814.697265625
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - example.sh
  - id: CreateReferenceFastaIndex
    in:
      - id: reference_fasta
        source: CreateReferenceFastaIndex.reference_fasta
      - id: monitoring_script
        source: CreateReferenceFastaIndex.monitoring_script
    out:
      - id: ref_index_output
      - id: monitoring_log
    run:
        class: CommandLineTool
        id: CreateReferenceFastaIndex
        inputs:
          - id: reference_fasta
            type: File
          - id: monitoring_script
            type:
              - File
              - 'null'
        outputs:
          - id: ref_index_output
            type: File
            outputBinding:
                glob: $(inputs.reference_fasta.basename + ".fai")
          - id: monitoring_log
            type: File
            outputBinding:
                glob: monitoring.log
        requirements:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/samtools:1.9
          - class: InitialWorkDirRequirement
            listing:
              - entryname: example.sh
                entry: |4

                    set -euo pipefail

                    # if the WDL/task contains a monitoring script as input
                    if [ ! -z "$(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)" ]; then
                      chmod a+x $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)
                      $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path) > monitoring.log &
                    else
                      echo "No monitoring script given as input" > monitoring.log &
                    fi

                    declare -a FASTA_TEMP=\$(mktemp fasta_XXXXXX)
                    cp $(inputs.reference_fasta.path) $FASTA_TEMP

                    # create a reference index
                    samtools faidx \
                        $FASTA_TEMP

                    mv $FASTA_TEMP.fai $(inputs.reference_fasta.basename + ".fai")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3337.860107421875
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - example.sh
