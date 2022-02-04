cwlVersion: v1.2
$graph:
  - class: CommandLineTool
    id: alignment_metrics
    inputs:
      - id: aligned_bam
        type: File
      - id: ref_fasta
        type: File
      - id: primers_bed
        type:
          - File
          - 'null'
      - id: machine_mem_gb
        type:
          - int
          - 'null'
      - id: docker
        default: quay.io/broadinstitute/viral-core:2.1.33
        type: string
    outputs:
      - id: wgs_metrics
        type: File
        outputBinding:
            glob: $(inputs.aligned_bam.basename.replace(/\.bam$/, '')  + '.raw_wgs_metrics.txt')
      - id: alignment_metrics
        type: File
        outputBinding:
            glob: $(inputs.aligned_bam.basename.replace(/\.bam$/, '')  + '.alignment_metrics.txt')
      - id: insert_size_metrics
        type: File
        outputBinding:
            glob: $(inputs.aligned_bam.basename.replace(/\.bam$/, '')  + '.insert_size_metrics.txt')
    requirements:
      - class: DockerRequirement
        dockerPull: quay.io/broadinstitute/viral-core:2.1.33
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                MEM_MB=\$(free -m | head -2 | tail -1 | awk '{print $4}')
                XMX=\$(echo "-Xmx"$MEM_MB"m")
                echo "Requesting $MEM_MB MB of RAM for Java"

                # requisite Picard fasta indexing
                cp "$(inputs.ref_fasta.path)" reference.fasta
                picard $XMX CreateSequenceDictionary -R reference.fasta

                # get Picard metrics and clean up the junky outputs
                picard $XMX CollectRawWgsMetrics \
                  -R reference.fasta \
                  -I "$(inputs.aligned_bam.path)" \
                  -O picard_raw.raw_wgs_metrics.txt
                grep -v \# picard_raw.raw_wgs_metrics.txt | grep . | head -2 > picard_clean.raw_wgs_metrics.txt

                picard $XMX CollectAlignmentSummaryMetrics \
                  -R reference.fasta \
                  -I "$(inputs.aligned_bam.path)" \
                  -O picard_raw.alignment_metrics.txt
                grep -v \# picard_raw.alignment_metrics.txt | grep . | head -4 > picard_clean.alignment_metrics.txt 

                picard $XMX CollectInsertSizeMetrics \
                  -I "$(inputs.aligned_bam.path)" \
                  -O picard_raw.insert_size_metrics.txt \
                  -H picard_raw.insert_size_metrics.pdf \
                  --INCLUDE_DUPLICATES true
                grep -v \# picard_raw.insert_size_metrics.txt | grep . | head -2 > picard_clean.insert_size_metrics.txt

                # prepend the sample name in order to facilitate tsv joining later
                SAMPLE=\$(samtools view -H "$(inputs.aligned_bam.path)" | grep ^@RG | perl -lape 's/^@RG.*SM:(\S+).*$/$1/' | sort | uniq)
                echo -e "sample_sanitized\tbam" > prepend.txt
                echo -e "$SAMPLE\t$(inputs.aligned_bam.basename.replace(/\.bam$/, '') )" >> prepend.txt
                paste prepend.txt picard_clean.raw_wgs_metrics.txt > "$(inputs.aligned_bam.basename.replace(/\.bam$/, '') )".raw_wgs_metrics.txt
                echo -e "$SAMPLE\t$(inputs.aligned_bam.basename.replace(/\.bam$/, '') )" >> prepend.txt
                echo -e "$SAMPLE\t$(inputs.aligned_bam.basename.replace(/\.bam$/, '') )" >> prepend.txt
                paste prepend.txt picard_clean.alignment_metrics.txt > "$(inputs.aligned_bam.basename.replace(/\.bam$/, '') )".alignment_metrics.txt
                echo -e "sample_sanitized\tbam" > prepend.txt
                echo -e "$SAMPLE\t$(inputs.aligned_bam.basename.replace(/\.bam$/, '') )" >> prepend.txt
                paste prepend.txt picard_clean.insert_size_metrics.txt > "$(inputs.aligned_bam.basename.replace(/\.bam$/, '') )".insert_size_metrics.txt

                # actually don't know how to do CollectTargetedPcrMetrics yet
                if [ -n "$(inputs.primers_bed === null ? "" : inputs.primers_bed.path)" ]; then
                  picard $XMX BedToIntervalList \
                    -I "$(inputs.primers_bed === null ? "" : inputs.primers_bed.path)" \
                    -O primers.interval.list \
                    -SD reference.dict
                fi
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 2
        ramMin: |-
            ${
            var unit = "GB";
            var value = parseInt(`${[inputs.machine_mem_gb, 13].find(element => element !== null) }`.match(/[0-9]+/g));
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
        outdirMin: 153600
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: plot_coverage
    inputs:
      - id: aligned_reads_bam
        type: File
      - id: sample_name
        type: string
      - id: skip_mark_dupes
        default: false
        type: boolean
      - id: plot_only_non_duplicates
        default: false
        type: boolean
      - id: bin_large_plots
        default: false
        type: boolean
      - id: binning_summary_statistic
        default: max
        type:
          - string
          - 'null'
      - id: plot_width_pixels
        default: 1100
        type:
          - int
          - 'null'
      - id: plot_height_pixels
        default: 850
        type:
          - int
          - 'null'
      - id: plot_pixels_per_inch
        default: 100
        type:
          - int
          - 'null'
      - id: max_coverage_depth
        type:
          - int
          - 'null'
      - id: base_q_threshold
        type:
          - int
          - 'null'
      - id: mapping_q_threshold
        type:
          - int
          - 'null'
      - id: read_length_threshold
        type:
          - int
          - 'null'
      - id: plotXLimits
        type:
          - string
          - 'null'
      - id: plotYLimits
        type:
          - string
          - 'null'
      - id: docker
        default: quay.io/broadinstitute/viral-core:2.1.33
        type: string
    outputs:
      - id: coverage_plot
        type: File
        outputBinding:
            glob: $(inputs.sample_name + '.coverage_plot.pdf')
      - id: coverage_tsv
        type: File
        outputBinding:
            glob: $(inputs.sample_name + '.coverage_plot.txt')
      - id: assembly_length
        type: int
        outputBinding:
            glob: $("assembly_length")
      - id: reads_aligned
        type: int
        outputBinding:
            glob: $("reads_aligned")
      - id: read_pairs_aligned
        type: int
        outputBinding:
            glob: $("read_pairs_aligned")
      - id: bases_aligned
        type: float
        outputBinding:
            loadContents: true
            glob: bases_aligned
            outputEval: $(parseFloat(self[0].contents))
      - id: mean_coverage
        type: float
        outputBinding:
            loadContents: true
            glob: mean_coverage
            outputEval: $(parseFloat(self[0].contents))
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

                samtools view -c $(inputs.aligned_reads_bam.path) | tee reads_aligned
                if [ "\$(cat reads_aligned)" != "0" ]; then
                  samtools index -@ "\$(nproc)" "$(inputs.aligned_reads_bam.path)"

                  PLOT_DUPE_OPTION=""
                  if [[ "$(inputs.skip_mark_dupes)" != "true" ]]; then
                    PLOT_DUPE_OPTION="$(inputs.plot_only_non_duplicates ? "--plotOnlyNonDuplicates" : "")"
                  fi

                  BINNING_OPTION="$(inputs.bin_large_plots ? "--binLargePlots" : "")"

                  # plot coverage
                  reports.py plot_coverage \
                    "$(inputs.aligned_reads_bam.path)" \
                    "$(inputs.sample_name).coverage_plot.pdf" \
                    --outSummary "$(inputs.sample_name).coverage_plot.txt" \
                    --plotFormat pdf \
                    --plotWidth 1100 \
                    --plotHeight 850 \
                    --plotDPI 100 \
                    $(inputs.max_coverage_depth === null ? "" : "-m " + inputs.max_coverage_depth) \
                    $(inputs.base_q_threshold === null ? "" : "-q " + inputs.base_q_threshold) \
                    $(inputs.mapping_q_threshold === null ? "" : "-Q " + inputs.mapping_q_threshold) \
                    $(inputs.read_length_threshold === null ? "" : "-l " + inputs.read_length_threshold) \
                    $(inputs.plotXLimits === null ? "" : "--plotXLimits " + inputs.plotXLimits) \
                    $(inputs.plotYLimits === null ? "" : "--plotYLimits " + inputs.plotYLimits) \
                    $PLOT_DUPE_OPTION \
                    $BINNING_OPTION \
                    --binningSummaryStatistic $(inputs.binning_summary_statistic) \
                    --plotTitle "$(inputs.sample_name) coverage plot" \
                    --loglevel=DEBUG

                else
                  touch $(inputs.sample_name).coverage_plot.pdf $(inputs.sample_name).coverage_plot.txt
                fi

                # collect figures of merit
                set +o pipefail # grep will exit 1 if it fails to find the pattern
                samtools view -H $(inputs.aligned_reads_bam.path) | perl -n -e'/^@SQ.*LN:(\d+)/ && print "$1\n"' |  python -c "import sys; print(sum(int(x) for x in sys.stdin))" | tee assembly_length
                # report only primary alignments 260=exclude unaligned reads and secondary mappings
                samtools view -h -F 260 $(inputs.aligned_reads_bam.path) | samtools flagstat - | tee $(inputs.sample_name).flagstat.txt
                grep properly $(inputs.sample_name).flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
                samtools view $(inputs.aligned_reads_bam.path) | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
                python -c "print (float("\$(cat bases_aligned)")/"\$(cat assembly_length)") if "\$(cat assembly_length)">0 else print(0)" > mean_coverage
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 2
        ramMin: 6675.72021484375
        outdirMin: 384000
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: coverage_report
    inputs:
      - id: mapped_bams
        type:
            items: File
            type: array
      - id: mapped_bam_idx
        type:
            items: File
            type: array
      - id: out_report_name
        default: coverage_report.txt
        type: string
      - id: docker
        default: quay.io/broadinstitute/viral-core:2.1.33
        type: string
    outputs:
      - id: coverage_report
        type: File
        outputBinding:
            glob: $(inputs.out_report_name)
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

                reports.py --version | tee VERSION
                reports.py coverage_only \
                  $(inputs.mapped_bams.map(function(el) {return el.path}).join(" ")) \
                  $(inputs.out_report_name) \
                  --loglevel DEBUG
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 2
        ramMin: 1907.3486328125
        outdirMin: 384000
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: assembly_bases
    inputs:
      - id: fasta
        type: File
      - id: docker
        default: ubuntu
        type: string
    outputs:
      - id: assembly_length
        type: int
        outputBinding:
            glob: $("assembly_length")
      - id: assembly_length_unambiguous
        type: int
        outputBinding:
            glob: $("assembly_length_unambiguous")
    requirements:
      - class: DockerRequirement
        dockerPull: ubuntu
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                grep -v '^>' "$(inputs.fasta.path)" | tr -d '\n' | wc -c | tee assembly_length
                grep -v '^>' "$(inputs.fasta.path)" | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 1
        ramMin: 953.67431640625
        outdirMin: 51200
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: fastqc
    inputs:
      - id: reads_bam
        type: File
      - id: docker
        default: quay.io/broadinstitute/viral-core:2.1.33
        type: string
    outputs:
      - id: fastqc_html
        type: File
        outputBinding:
            glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '_fastqc.html')
      - id: fastqc_zip
        type: File
        outputBinding:
            glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '_fastqc.zip')
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
                reports.py --version | tee VERSION
                reports.py fastqc $(inputs.reads_bam.path) $(inputs.reads_bam.basename.replace(/\.bam$/, '') )_fastqc.html --out_zip $(inputs.reads_bam.basename.replace(/\.bam$/, '') )_fastqc.zip
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 1
        ramMin: 1907.3486328125
        outdirMin: 384000
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
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
            glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '.count.' +
                inputs.ref_db.basename.replace(/\.fasta$/, '')  + '.txt')
      - id: report_top_hits
        type: File
        outputBinding:
            glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '.count.' +
                inputs.ref_db.basename.replace(/\.fasta$/, '')  + '.top_' + inputs.topNHits
                + '_hits.txt')
      - id: top_hit_id
        type: string
        outputBinding:
            loadContents: true
            glob: $(inputs.reads_bam.basename.replace(/\.bam$/, '')  + '.count.' +
                inputs.ref_db.basename.replace(/\.fasta$/, '')  + '.top.txt')
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
        ramMin: 14305.11474609375
        outdirMin: 384000
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
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
      - class: DockerRequirement
        dockerPull: quay.io/broadinstitute/viral-core:2.1.33
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
      - class: ResourceRequirement
        coresMin: 8
        ramMin: 6675.72021484375
        outdirMin: 102400
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: aggregate_metagenomics_reports
    inputs:
      - id: kraken_summary_reports
        type:
            items: File
            type: array
      - id: aggregate_taxon_heading_space_separated
        default: Viruses
        type: string
      - id: aggregate_taxlevel_focus
        default: species
        type: string
      - id: aggregate_top_N_hits
        default: 5
        type: int
      - id: docker
        default: quay.io/broadinstitute/viral-classify:2.1.16.0
        type: string
    outputs:
      - id: krakenuniq_aggregate_taxlevel_summary
        type: File
        outputBinding:
            glob: $('aggregate_taxa_summary_' + inputs.aggregate_taxon_heading_space_separated.replace("
                ", "_")  + '_by_' + inputs.aggregate_taxlevel_focus + '_top_' + inputs.aggregate_top_N_hits
                + '_by_sample.csv')
      - id: viralngs_version
        type: string
        outputBinding:
            loadContents: true
            glob: VERSION
            outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
    requirements:
      - class: DockerRequirement
        dockerPull: quay.io/broadinstitute/viral-classify:2.1.16.0
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -ex -o pipefail

                metagenomics.py --version | tee VERSION
                metagenomics.py taxlevel_summary \
                  $(inputs.kraken_summary_reports.map(function(el) {return el.path}).join(" ")) \
                  --csvOut aggregate_taxa_summary_$(inputs.aggregate_taxon_heading_space_separated.replace(" ", "_") )_by_$(inputs.aggregate_taxlevel_focus)_top_$(inputs.aggregate_top_N_hits)_by_sample.csv \
                  --noHist \
                  --taxHeading $(inputs.aggregate_taxon_heading_space_separated) \
                  --taxlevelFocus $(inputs.aggregate_taxlevel_focus) \
                  --zeroFill --includeRoot --topN $(inputs.aggregate_top_N_hits) \
                  --loglevel=DEBUG
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 1
        ramMin: 2861.02294921875
        outdirMin: 51200
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: MultiQC
    inputs:
      - id: input_files
        default: []
        type:
            items: File
            type: array
      - id: force
        default: false
        type: boolean
      - id: full_names
        default: false
        type: boolean
      - id: title
        type:
          - string
          - 'null'
      - id: comment
        type:
          - string
          - 'null'
      - id: file_name
        type:
          - string
          - 'null'
      - id: out_dir
        default: ./multiqc-output
        type: string
      - id: template
        type:
          - string
          - 'null'
      - id: tag
        type:
          - string
          - 'null'
      - id: ignore_analysis_files
        type:
          - string
          - 'null'
      - id: ignore_sample_names
        type:
          - string
          - 'null'
      - id: sample_names
        type:
          - File
          - 'null'
      - id: exclude_modules
        type:
          - items: string
            type: array
          - 'null'
      - id: module_to_use
        type:
          - items: string
            type: array
          - 'null'
      - id: data_dir
        default: false
        type: boolean
      - id: no_data_dir
        default: false
        type: boolean
      - id: output_data_format
        type:
          - string
          - 'null'
      - id: zip_data_dir
        default: false
        type: boolean
      - id: export
        default: false
        type: boolean
      - id: flat
        default: false
        type: boolean
      - id: interactive
        default: true
        type: boolean
      - id: lint
        default: false
        type: boolean
      - id: pdf
        default: false
        type: boolean
      - id: megaQC_upload
        default: false
        type: boolean
      - id: config
        type:
          - File
          - 'null'
      - id: config_yaml
        type:
          - string
          - 'null'
      - id: docker
        default: quay.io/biocontainers/multiqc:1.8--py_2
        type: string
    outputs:
      - id: multiqc_report
        type: File
        outputBinding:
            glob: "$(inputs.out_dir + '/' + inputs.file_name ? [inputs.file_name].find(element\
                \ => element !== null) .split('/').reverse()[0].replace(/\\.html$/,\
                \ '') : \"multiqc\" + '.html')"
      - id: multiqc_data_dir_tarball
        type: File
        outputBinding:
            glob: "$(inputs.file_name ? [inputs.file_name].find(element => element\
                \ !== null) .split('/').reverse()[0].replace(/\\.html$/, '') : \"\
                multiqc\" + '_data.tar.gz')"
    requirements:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/multiqc:1.8--py_2
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                      set -ex -o pipefail

                      echo "$(inputs.input_files.map(function(el) {return el.path}).join("
                "))" > input-filenames.txt
                      echo "" >> input-filenames.txt

                      multiqc \
                      --file-list input-filenames.txt \
                      --dirs \
                      --outdir "$(inputs.out_dir)" \
                      $(inputs.force ? "--force" : "") \
                      $(inputs.full_names ? "--fullnames" : "") \
                      $(inputs.title === null ? "" : "--title " + inputs.title) \
                      $(inputs.comment === null ? "" : "--comment " + inputs.comment) \
                      $(inputs.file_name === null ? "" : "--filename " + inputs.file_name) \
                      $(inputs.template === null ? "" : "--template " + inputs.template) \
                      $(inputs.tag === null ? "" : "--tag " + inputs.tag) \
                      $(inputs.ignore_analysis_files === null ? "" : "--ignore " + inputs.ignore_analysis_files) \
                      $(inputs.ignore_sample_names === null ? "" : "--ignore-samples" + inputs.ignore_sample_names) \
                      $(inputs.sample_names === null ? "" : "--sample-names " + inputs.sample_names.path) \
                      $(inputs.exclude_modules === null ? "" : "--exclude ")$([inputs.exclude_modules, []].find(element => element !== null) .join(" --exclude ")) \
                      $(inputs.module_to_use === null ? "" : "--module ")$([inputs.module_to_use, []].find(element => element !== null) .join(" --module ")) \
                      $(inputs.data_dir ? "--data-dir" : "") \
                      $(inputs.no_data_dir ? "--no-data-dir" : "") \
                      $(inputs.output_data_format === null ? "" : "--data-format " + inputs.output_data_format) \
                      $(inputs.zip_data_dir ? "--zip-data-dir" : "") \
                      $(inputs["export"] ? "--export" : "") \
                      $(inputs.flat ? "--flat" : "") \
                      $(inputs.interactive ? "--interactive" : "") \
                      $(inputs.lint ? "--lint" : "") \
                      $(inputs.pdf ? "--pdf" : "") \
                      $(inputs.megaQC_upload ? "" : "--no-megaqc-upload") \
                      $(inputs.config === null ? "" : "--config " + inputs.config.path) \
                      $(inputs.config_yaml === null ? "" : "--cl-config " + inputs.config_yaml)

                      if [ -z "$(inputs.file_name)" ]; then
                        mv "$(inputs.out_dir)/$(inputs.file_name ? [inputs.file_name].find(element => element !== null) .split('/').reverse()[0].replace(/\.html$/, '') : "multiqc")_report.html" "$(inputs.out_dir)/$(inputs.file_name ? [inputs.file_name].find(element => element !== null) .split('/').reverse()[0].replace(/\.html$/, '') : "multiqc").html"
                      fi

                      tar -c "$(inputs.out_dir)/$(inputs.file_name ? [inputs.file_name].find(element => element !== null) .split('/').reverse()[0].replace(/\.html$/, '') : "multiqc")_data" | gzip -c > "$(inputs.file_name ? [inputs.file_name].find(element => element !== null) .split('/').reverse()[0].replace(/\.html$/, '') : "multiqc")_data.tar.gz"
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 16
        ramMin: 7629.39453125
        outdirMin: 384000
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: compare_two_genomes
    inputs:
      - id: genome_one
        type: File
      - id: genome_two
        type: File
      - id: out_basename
        type: string
      - id: docker
        default: quay.io/broadinstitute/viral-assemble:2.1.16.1
        type: string
    outputs:
      - id: comparison_table
        type: File
        outputBinding:
            glob: $(inputs.out_basename + '.txt')
      - id: max_ram_gb
        type: int
        outputBinding:
            glob: $(Math.ceil("MEM_BYTES"/1000000000) )
      - id: runtime_sec
        type: int
        outputBinding:
            glob: $(Math.ceil("UPTIME_SEC") )
      - id: cpu_load
        type: string
        outputBinding:
            loadContents: true
            glob: CPU_LOAD
            outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
      - id: viralngs_version
        type: string
        outputBinding:
            loadContents: true
            glob: VERSION
            outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
    requirements:
      - class: DockerRequirement
        dockerPull: quay.io/broadinstitute/viral-assemble:2.1.16.1
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -ex -o pipefail
                assembly.py --version | tee VERSION
                assembly.py alignment_summary "$(inputs.genome_one.path)" "$(inputs.genome_two.path)" --outfileName "$(inputs.out_basename).txt" --printCounts --loglevel=DEBUG
                cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
                cat /proc/loadavg > CPU_LOAD
                cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
      - class: ResourceRequirement
        coresMin: 2
        ramMin: 2861.02294921875
        outdirMin: 51200
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
