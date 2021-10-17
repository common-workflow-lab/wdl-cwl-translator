class: CommandLineTool
id: ValidateMetrics
inputs:
  - id: cell_metrics
    type: File
  - id: gene_metrics
    type: File
  - id: expected_cell_metric_hash
    type: string
  - id: expected_gene_metric_hash
    type: string
outputs:
  - id: result
    type: string
    outputBinding:
        loadContents: true
        glob: result.txt
        outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
requirements:
  - class: DockerRequirement
    dockerPull: ubuntu:18.04
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -eo pipefail

            cacheInvalidationRandomString=4

            # check matrix row and column indexes files hash
            gene_metric_hash=\$(zcat "$(inputs.gene_metrics.path)" | md5sum | awk '{print $1}')
            cell_metric_hash=\$(zcat "$(inputs.cell_metrics.path)" | md5sum | awk '{print $1}')

            fail=false

            if [ "$gene_metric_hash" == "$(inputs.expected_gene_metric_hash)" ]; then
                echo Computed and expected gene metrics match \( "$gene_metric_hash" \)
            else
                echo Computed \( "$gene_metric_hash" \) and expected \( "$(inputs.expected_gene_metric_hash)" \) gene checksums do not match
                fail=true
            fi

            if [ "$cell_metric_hash" == "$(inputs.expected_cell_metric_hash)" ]; then
                echo Computed and expected cell metrics match \( "$cell_metric_hash" \)
            else
                echo Computed \( "$cell_metric_hash" \) and expected \( "$(inputs.expected_cell_metric_hash)" \) cell metrics hashes do not match
                fail=true
            fi

            if [ $fail == "true" ]; then
                printf FAIL > result.txt
            else
                printf PASS > result.txt
            fi
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    ramMin: 0
  - class: ResourceRequirement
    coresMin: 1
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
