cwlVersion: v1.2
id: BuildCembaReferences
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
inputs:
  - id: reference_fasta
    type: File
  - id: monitoring_script
    type:
      - File
      - 'null'
steps:
  - id: Convert
    in:
      - id: fasta_input
        source: reference_fasta
      - id: monitoring_script
        source: monitoring_script
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
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    # if the WDL/task contains a monitoring script as input
                    if [ ! -z "$(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)" ]; then
                      chmod a+x $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path)
                      $(inputs.monitoring_script === null ? "" : inputs.monitoring_script.path) > monitoring.log &
                    else
                      echo "No monitoring script given as input" > monitoring.log &
                    fi

                    python3 /build_bisulfite_references.py \
                      --input-fasta $(inputs.fasta_input.path) \
                      --forward-convert-out $("genome_mfa.CT_conversion.fa") \
                      --reverse-convert-out $("genome_mfa.GA_conversion.fa")
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/bisulfite-references:1.0
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3337.860107421875
            outdirMin: '$((Math.ceil(3.5 * (function(size_of=0){inputs.fasta_input.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3 < 1 ? 1 : (function(size_of=0){inputs.fasta_input.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: IndexForward
    in:
      - id: fasta_input
        source: Convert/fwd_converted_reference_fasta_output
      - id: index_prefix
        default: BS_CT
      - id: monitoring_script
        source: monitoring_script
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
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
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
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/bowtie2:2.3.4.3
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 6675.72021484375
            outdirMin: '$((Math.ceil(3 * (function(size_of=0){inputs.fasta_input.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3 < 1 ? 1 : (function(size_of=0){inputs.fasta_input.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: IndexReverse
    in:
      - id: fasta_input
        source: Convert/rev_converted_reference_fasta_output
      - id: index_prefix
        default: BS_GA
      - id: monitoring_script
        source: monitoring_script
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
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
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
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/bowtie2:2.3.4.3
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 6675.72021484375
            outdirMin: '$((Math.ceil(3 * (function(size_of=0){inputs.fasta_input.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3 < 1 ? 1 : (function(size_of=0){inputs.fasta_input.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: CreateReferenceDictionary
    in:
      - id: reference_fasta
        source: reference_fasta
      - id: monitoring_script
        source: monitoring_script
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
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
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
                    java -Xmx3500m -jar /picard-tools/picard.jar CreateSequenceDictionary \
                      REFERENCE=$(inputs.reference_fasta.path) \
                      OUTPUT=$(inputs.reference_fasta.basename.replace(/\.fa$/, '') .split('/').reverse()[0].replace(/\.fasta$/, '') + ".dict")
                    sed -i "s=\$(dirname $(inputs.reference_fasta.path))/==g" $(inputs.reference_fasta.basename.replace(/\.fa$/, '') .split('/').reverse()[0].replace(/\.fasta$/, '') + ".dict")  # for reproducibility

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/picard:2.18.23
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 4000.0
            outdirMin: '$((Math.ceil(2 * (function(size_of=0){inputs.reference_fasta.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3 < 1 ? 1 : (function(size_of=0){inputs.reference_fasta.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: CreateReferenceFastaIndex
    in:
      - id: reference_fasta
        source: reference_fasta
      - id: monitoring_script
        source: monitoring_script
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
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
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
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/samtools:1.9
          - class: ResourceRequirement
            coresMin: 1
            ramMin: 3337.860107421875
            outdirMin: '$((Math.ceil(2.25 * (function(size_of=0){inputs.reference_fasta.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3 < 1 ? 1 : (function(size_of=0){inputs.reference_fasta.path.forEach(function(element){
                if (element) {size_of += element.size}})}) / 1000^3) ) * 1024)'
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
outputs:
  - id: BuildCembaReferences.reference_fasta_dict
    outputSource: CreateReferenceDictionary/ref_dict_output
    type: File
  - id: BuildCembaReferences.reference_fasta_index
    outputSource: CreateReferenceFastaIndex/ref_index_output
    type: File
  - id: BuildCembaReferences.fwd_converted_reference_fasta
    outputSource: Convert/fwd_converted_reference_fasta_output
    type: File
  - id: BuildCembaReferences.rev_converted_reference_fasta
    outputSource: Convert/rev_converted_reference_fasta_output
    type: File
  - id: BuildCembaReferences.fwd_bowtie2_index_files
    outputSource: IndexForward/bowtie2_index_files
    type:
        items: File
        type: array
  - id: BuildCembaReferences.rev_bowtie2_index_files
    outputSource: IndexReverse/bowtie2_index_files
    type:
        items: File
        type: array
