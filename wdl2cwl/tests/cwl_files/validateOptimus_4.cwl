cwlVersion: v1.2
id: ValidateMetrics
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
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
hints:
  - class: DockerRequirement
    dockerPull: ubuntu:18.04
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 953.67431640625
    outdirMin: $((Math.ceil((function(size_of=0){inputs.cell_metrics.forEach(function(element){
        if (element) {size_of += element.size}})}) / 1000^3 + (function(size_of=0){inputs.gene_metrics.forEach(function(element){
        if (element) {size_of += element.size}})}) / 1000^3 * 1.1) ) * 1024)
inputs:
  - id: cell_metrics
    type: File
  - id: gene_metrics
    type: File
  - id: expected_cell_metric_hash
    type: string
  - id: expected_gene_metric_hash
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: result
    type: string
    outputBinding:
        loadContents: true
        glob: result.txt
        outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
