cwlVersion: v1.2
id: align_and_count_multiple_report
class: Workflow
doc: Count the number of times reads map to provided reference sequences. Useful for
    counting spike-ins, etc.
requirements:
  - class: ScatterFeatureRequirement
inputs:
  - id: reads_unmapped_bams
    doc: Unaligned reads in BAM format
    type:
        items: File
        type: array
  - id: ref_db
    doc: File containing sequences against which reads should me aligned and counted
    type: File
  - id: align_and_count.topNHits
    default: 3
    type: int
  - id: align_and_count.machine_mem_gb
    type:
      - int
      - 'null'
  - id: align_and_count.docker
    default: quay.io/broadinstitute/viral-core:2.1.33
    type: string
  - id: align_and_count_summary.output_prefix
    default: count_summary
    type: string
  - id: align_and_count_summary.docker
    default: quay.io/broadinstitute/viral-core:2.1.33
    type: string
  - id: align_and_count_summary_top_hits.docker
    default: quay.io/broadinstitute/viral-core:2.1.33
    type: string
steps:
  - id: align_and_count
    in:
      - id: reads_bam
        source: reads_unmapped_bams
      - id: ref_db
        source: ref_db
      - id: topNHits
        source: align_and_count.topNHits
      - id: machine_mem_gb
        source: align_and_count.machine_mem_gb
      - id: docker
        source: align_and_count.docker
    out:
      - id: report
      - id: report_top_hits
      - id: top_hit_id
      - id: viralngs_version
    run:
        class: CommandLineTool
        id: align_and_count
        inputs:
          - id: reads_bam
            type: File
          - id: ref_db
            type: File
          - id: topNHits
            default: 3
            type: int
          - id: machine_mem_gb
            type:
              - int
              - 'null'
          - id: docker
            default: quay.io/broadinstitute/viral-core:2.1.33
            type: string
        outputs:
          - id: report
            type: File
            outputBinding:
                glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '.count.'
                    + inputs.ref_db.basename.replace(/\.fasta$/, '')  + '.txt')
          - id: report_top_hits
            type: File
            outputBinding:
                glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '.count.'
                    + inputs.ref_db.basename.replace(/\.fasta$/, '')  + '.top_' +
                    inputs.topNHits + '_hits.txt')
          - id: top_hit_id
            type: string
            outputBinding:
                loadContents: true
                glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '.count.'
                    + inputs.ref_db.basename.replace(/\.fasta$/, '')  + '.top.txt')
                outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
          - id: viralngs_version
            type: string
            outputBinding:
                loadContents: true
                glob: VERSION
                outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -ex -o pipefail

                    read_utils.py --version | tee VERSION

                    ln -s "$(inputs.reads_bam.path)" "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).bam"
                    read_utils.py minimap2_idxstats \
                      "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).bam" \
                      "$(inputs.ref_db.path)" \
                      --outStats "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).count.$(inputs.ref_db.basename.replace(/\.fasta$/, '') ).txt.unsorted" \
                      --loglevel=DEBUG

                    sort -b -r -n -k3 "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).count.$(inputs.ref_db.basename.replace(/\.fasta$/, '') ).txt.unsorted" > "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).count.$(inputs.ref_db.basename.replace(/\.fasta$/, '') ).txt"
                    head -n $(inputs.topNHits) "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).count.$(inputs.ref_db.basename.replace(/\.fasta$/, '') ).txt" > "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).count.$(inputs.ref_db.basename.replace(/\.fasta$/, '') ).top_$(inputs.topNHits)_hits.txt"
                    head -1 "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).count.$(inputs.ref_db.basename.replace(/\.fasta$/, '') ).txt" | cut -f 1 > "$(inputs.reads_bam.basename.replace(/\.bam$/, '') ).count.$(inputs.ref_db.basename.replace(/\.fasta$/, '') ).top.txt"
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/viral-core:2.1.33
          - class: ResourceRequirement
            coresMin: 4
            ramMin: |-
                ${
                var unit = "GB";
                var value = parseInt(`${[inputs.machine_mem_gb, 15].find(function(element) { return element !== null }) }`.match(/[0-9]+/g));
                var memory = "";
                if(unit==="KiB") memory = value/1024;
                else if(unit==="MiB") memory = value;
                else if(unit==="GiB") memory = value*1024;
                else if(unit==="TiB") memory = value*1024*1024;
                else if(unit==="B") memory = value/(1024*1024);
                else if(unit==="KB" || unit==="K") memory = (value*1000)/(1024*1024);
                else if(unit==="MB" || unit==="M") memory = (value*(1000*1000))/(1024*1024);
                else if(unit==="GB" || unit==="G") memory = (value*(1000*1000*1000))/(1024*1024);
                else if(unit==="TB" || unit==="T") memory = (value*(1000*1000*1000*1000))/(1024*1024);
                else throw "Unknown units: " + unit;
                return parseInt(memory);
                }
            outdirMin: 384000
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
    scatter: reads_bam
  - id: align_and_count_summary
    in:
      - id: counts_txt
        source: align_and_count/report
      - id: output_prefix
        source: align_and_count_summary.output_prefix
      - id: docker
        source: align_and_count_summary.docker
    out:
      - id: count_summary
      - id: viralngs_version
    run:
        class: CommandLineTool
        id: align_and_count_summary
        inputs:
          - id: counts_txt
            type:
                items: File
                type: array
          - id: output_prefix
            default: count_summary
            type: string
          - id: docker
            default: quay.io/broadinstitute/viral-core:2.1.33
            type: string
        outputs:
          - id: count_summary
            type: File
            outputBinding:
                glob: $(inputs.output_prefix + '.tsv')
          - id: viralngs_version
            type: string
            outputBinding:
                loadContents: true
                glob: VERSION
                outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -ex -o pipefail

                    reports.py --version | tee VERSION
                    reports.py aggregate_alignment_counts $(inputs.counts_txt.map(function(el) {return el.path}).join(" ")) "$(inputs.output_prefix)".tsv --loglevel=DEBUG
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/viral-core:2.1.33
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 6675.72021484375
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
  - id: align_and_count_summary_top_hits
    in:
      - id: counts_txt
        source: align_and_count/report_top_hits
      - id: output_prefix
        default: count_summary_top_hits
      - id: docker
        source: align_and_count_summary_top_hits.docker
    out:
      - id: count_summary
      - id: viralngs_version
    run:
        class: CommandLineTool
        id: align_and_count_summary
        inputs:
          - id: counts_txt
            type:
                items: File
                type: array
          - id: output_prefix
            default: count_summary
            type: string
          - id: docker
            default: quay.io/broadinstitute/viral-core:2.1.33
            type: string
        outputs:
          - id: count_summary
            type: File
            outputBinding:
                glob: $(inputs.output_prefix + '.tsv')
          - id: viralngs_version
            type: string
            outputBinding:
                loadContents: true
                glob: VERSION
                outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -ex -o pipefail

                    reports.py --version | tee VERSION
                    reports.py aggregate_alignment_counts $(inputs.counts_txt.map(function(el) {return el.path}).join(" ")) "$(inputs.output_prefix)".tsv --loglevel=DEBUG
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/viral-core:2.1.33
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 6675.72021484375
            outdirMin: 102400
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
outputs:
  - id: align_and_count_multiple_report.report
    outputSource: align_and_count_summary/count_summary
    type: File
  - id: align_and_count_multiple_report.report_top_hits
    outputSource: align_and_count_summary_top_hits/count_summary
    type: File
  - id: align_and_count_multiple_report.viral_core_version
    outputSource: align_and_count_summary/viralngs_version
    type: string
