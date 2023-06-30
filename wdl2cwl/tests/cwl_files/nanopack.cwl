cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: NanoPlot
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputDir + inputs.outputPrefix))"
                NanoPlot \
                --threads $(inputs.threads) \
                --outdir $(inputs.outputDir) \
                --prefix $(inputs.outputPrefix) \
                $(inputs.outputTsvStats ? "--tsv_stats" : "") \
                $(inputs.dropOutliers ? "--drop_outliers" : "") \
                $(inputs.logLengths ? "--loglength" : "") \
                --format $(inputs.format) \
                $(inputs.showN50 ? "--N50" : "--no-N50") \
                $(inputs.maxLength === null ? "" : "--maxlength " + inputs.maxLength) \
                $(inputs.minLength === null ? "" : "--minlength " + inputs.minLength) \
                $(inputs.minQual === null ? "" : "--minqual " + inputs.minQual) \
                $(inputs.readType === null ? "" : "--readtype " + inputs.readType) \
                $({ "fastq": "--fastq ", "fasta": "--fasta ", "fastq_rich": "--fastq_rich ", "fastq_minimal": "--fastq_minimal ", "summary": "--summary ", "bam": "--bam ", "ubam": "--ubam ", "cram": "--cram ", "pickle": "--pickle ", "feather": "--feather " }[inputs.inputFileType] + inputs.inputFile.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0
      - class: ResourceRequirement
        coresMin: $(inputs.threads)
        ramMin: |-
            ${
            var unit = inputs.memory.match(/[a-zA-Z]+/g).join("");
            var value = parseInt(`${inputs.memory}`.match(/[0-9]+/g));
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
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: inputFile
        doc: The input file.
        type: File
      - id: inputFileType
        doc: The format of the read file.
        type: string
      - id: outputDir
        doc: Output directory path.
        type: string
      - id: outputPrefix
        doc: Output file prefix.
        type: string
      - id: outputPath
        doc: Combination of the outputDir & outputPrefix strings.
        type:
          - string
          - 'null'
      - id: outputTsvStats
        doc: Output the stats file as a properly formatted TSV.
        default: true
        type: boolean
      - id: dropOutliers
        doc: Drop outlier reads with extreme long length.
        default: false
        type: boolean
      - id: logLengths
        doc: Additionally show logarithmic scaling of lengths in plots.
        default: false
        type: boolean
      - id: format
        doc: Specify the output format of the plots.
        default: png
        type: string
      - id: showN50
        doc: Show the N50 mark in the read length histogram.
        default: true
        type: boolean
      - id: title
        doc: Add a title to all plots, requires quoting if using spaces.
        type:
          - string
          - 'null'
      - id: maxLength
        doc: Hide reads longer than length specified.
        type:
          - int
          - 'null'
      - id: minLength
        doc: Hide reads shorter than length specified.
        type:
          - int
          - 'null'
      - id: minQual
        doc: Drop reads with an average quality lower than specified.
        type:
          - int
          - 'null'
      - id: readType
        doc: Which read type to extract information about from summary. Options are
            1D, 2D, 1D2
        type:
          - string
          - 'null'
      - id: threads
        doc: The number of threads to be used.
        default: 2
        type: int
      - id: memory
        doc: The amount of memory available to the job.
        default: 2G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 15
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: dynamicHistogram
        doc: Dynamic histogram of read length.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "Dynamic_Histogram_Read_length.html")
      - id: readLengthHistogram
        doc: Histogram of read length.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "HistogramReadlength.png")
      - id: logScaleReadLengthHistogram
        doc: Histogram of read lengths after log transformation.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "LogTransformed_HistogramReadlength.png")
      - id: report
        doc: Html summary report.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "NanoPlot-report.html")
      - id: weightedHistogram
        doc: Weighted histogram of read lengths.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "Weighted_HistogramReadlength.png")
      - id: weightedLogScaleHistogram
        doc: Weighted histogram of read lengths after log transformation.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "Weighted_LogTransformed_HistogramReadlength.png")
      - id: yieldByLength
        doc: Cumulative yield plot.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "Yield_By_Length.png")
      - id: lengthVsQualityScatterPlotDot
        doc: Read lengths vs average read quality plot.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "LengthvsQualityScatterPlot_dot.png")
      - id: lengthVsQualityScatterPlotKde
        doc: Read lengths vs average read quality plot.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "LengthvsQualityScatterPlot_kde.png")
      - id: stats
        doc: NanoStats report.
        type:
          - File
          - 'null'
        outputBinding:
            glob: $(inputs.outputDir + inputs.outputPrefix + "NanoStats.txt")
  - cwlVersion: v1.2
    id: NanoQc
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputDir))"
                nanoQC \
                --outdir $(inputs.outputDir) \
                $(inputs.directRna ? "--rna" : "") \
                $(inputs.minLength === null ? "" : "--minlen " + inputs.minLength) \
                $(inputs.inputFile.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: quay.io/biocontainers/nanoqc:0.9.4--py_0
      - class: ResourceRequirement
        ramMin: |-
            ${
            var unit = inputs.memory.match(/[a-zA-Z]+/g).join("");
            var value = parseInt(`${inputs.memory}`.match(/[0-9]+/g));
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
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    inputs:
      - id: inputFile
        doc: The input file.
        type: File
      - id: outputDir
        doc: Output directory path.
        type: string
      - id: directRna
        doc: Fastq is from direct RNA-seq and contains U nucleotides.
        default: false
        type: boolean
      - id: minLength
        doc: Filters the reads on a minimal length of the given range. Also plots
            the given length/2 of the begin and end of the reads.
        type:
          - int
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 2G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 15
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: quay.io/biocontainers/nanoqc:0.9.4--py_0
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: report
        doc: Html summary report.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "nanoQC.html")
      - id: log
        doc: Progress report.
        type: File
        outputBinding:
            glob: $(inputs.outputDir + "NanoQC.log")
