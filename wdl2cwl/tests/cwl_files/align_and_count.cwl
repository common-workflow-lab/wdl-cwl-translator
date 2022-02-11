class: Workflow
id: align_and_count_report
doc: Align reads to reference with minimap2 and count the number of hits. Results
    are returned in the format of 'samtools idxstats'.
inputs:
  - id: align_and_count.reads_bam
    type: File
  - id: align_and_count.ref_db
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
outputs:
  - id: align_and_count_report.report
    outputSource: reports.align_and_count/report
    type: File
  - id: align_and_count_report.report_top_hits
    outputSource: reports.align_and_count/report_top_hits
    type: File
  - id: align_and_count_report.viral_core_version
    outputSource: reports.align_and_count/viralngs_version
    type: string
cwlVersion: v1.2
steps:
  - id: reports.align_and_count
    in:
      - id: reads_bam
        source: align_and_count.reads_bam
      - id: ref_db
        source: align_and_count.ref_db
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
          - class: DockerRequirement
            dockerPull: quay.io/broadinstitute/viral-core:2.1.33
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
                return parseInt(memory);
                }
            outdirMin: 384000
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
