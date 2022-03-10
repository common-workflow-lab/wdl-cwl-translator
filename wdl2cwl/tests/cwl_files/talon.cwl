cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: CreateAbundanceFileFromDatabase
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_abundance \
                --db=$(inputs.databaseFile.path) \
                -a $(inputs.annotationVersion) \
                -b $(inputs.genomeBuild) \
                --o=$(inputs.outputPrefix) \
                $(inputs.whitelistFile === null ? "" : "--whitelist=" + inputs.whitelistFile.path) \
                $(inputs.datasetsFile === null ? "" : "-d " + inputs.datasetsFile.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: databaseFile
        doc: Talon database.
        type: File
      - id: annotationVersion
        doc: Which annotation version to use.
        type: string
      - id: genomeBuild
        doc: Genome build to use.
        type: string
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: whitelistFile
        doc: Whitelist file of transcripts to include in the output.
        type:
          - File
          - 'null'
      - id: datasetsFile
        doc: A file indicating which datasets should be included.
        type:
          - File
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: abundanceFile
        doc: Abundance for each transcript in the talon database across datasets.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_talon_abundance.tsv")
  - cwlVersion: v1.2
    id: CreateGtfFromDatabase
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_create_GTF \
                --db=$(inputs.databaseFile.path) \
                -b $(inputs.genomeBuild) \
                -a $(inputs.annotationVersion) \
                --o=$(inputs.outputPrefix) \
                $(inputs.observedInDataset ? "--observed" : "") \
                $(inputs.whitelistFile === null ? "" : "--whitelist=" + inputs.whitelistFile.path) \
                $(inputs.datasetFile === null ? "" : "-d " + inputs.datasetFile.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: databaseFile
        doc: Talon database.
        type: File
      - id: genomeBuild
        doc: Genome build to use.
        type: string
      - id: annotationVersion
        doc: Which annotation version to use.
        type: string
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: observedInDataset
        doc: The output will only include transcripts that were observed at least
            once.
        default: false
        type: boolean
      - id: whitelistFile
        doc: Whitelist file of transcripts to include in the output.
        type:
          - File
          - 'null'
      - id: datasetFile
        doc: A file indicating which datasets should be included.
        type:
          - File
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: gtfFile
        doc: The genes, transcripts, and exons stored a talon database in gtf format.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_talon.gtf")
  - cwlVersion: v1.2
    id: FilterTalonTranscripts
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_filter_transcripts \
                --db=$(inputs.databaseFile.path) \
                -a $(inputs.annotationVersion) \
                $("--o=" + inputs.outputPrefix + "_whitelist.csv") \
                --maxFracA=$(inputs.maxFracA) \
                --minCount=$(inputs.minCount) \
                $(inputs.allowGenomic ? "--allowGenomic" : "") \
                --datasets=$(inputs.datasetsFile === null ? "" : inputs.datasetsFile.path) \
                --minDatasets=$(inputs.minDatasets)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: databaseFile
        doc: Talon database.
        type: File
      - id: annotationVersion
        doc: Which annotation version to use.
        type: string
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: maxFracA
        doc: Maximum fraction of As to allow in the window located immediately after
            any read assigned to a novel transcript.
        default: 0.5
        type: float
      - id: minCount
        doc: Number of minimum occurrences required for a novel transcript per dataset.
        default: 5
        type: int
      - id: allowGenomic
        doc: If this option is set, transcripts from the Genomic novelty category
            will be permitted in the output.
        default: false
        type: boolean
      - id: datasetsFile
        doc: Datasets to include.
        type:
          - File
          - 'null'
      - id: minDatasets
        doc: Minimum number of datasets novel transcripts must be found in.
        type:
          - int
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: transcriptWhitelist
        doc: Transcript whitelist produced from the talon database.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_whitelist.csv")
  - cwlVersion: v1.2
    id: GetReadAnnotations
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_fetch_reads \
                --db $(inputs.databaseFile.path) \
                --build $(inputs.genomeBuild) \
                --o $(inputs.outputPrefix) \
                $(inputs.datasetFile === null ? "" : "--datasets " + inputs.datasetFile.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: databaseFile
        doc: Talon database.
        type: File
      - id: genomeBuild
        doc: Genome build to use.
        type: string
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: datasetFile
        doc: A file indicating which datasets should be included.
        type:
          - File
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: readAnnotations
        doc: Read-specific annotation information from a talon database.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_talon_read_annot.tsv")
  - cwlVersion: v1.2
    id: GetSpliceJunctions
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_get_sjs \
                $({ "db": "--db", "gtf": "--gtf" }[inputs.inputFileType] + inputs.sjInformationFile.path) \
                --ref $(inputs.referenceGtf.path) \
                --mode $(inputs.runMode) \
                --outprefix $(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: sjInformationFile
        doc: Talon gtf file or database from which to extract exons/introns.
        type: File
      - id: inputFileType
        doc: The file type of sjInformationFile.
        default: db
        type: string
      - id: referenceGtf
        doc: Gtf reference file (ie gencode).
        type: File
      - id: runMode
        doc: Determines whether to include introns or exons in the output.
        default: intron
        type: string
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: spliceJunctions
        doc: File containing locations, novelty and transcript assignments of exons/introns.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_" + inputs.runMode + "s.tsv")
  - cwlVersion: v1.2
    id: InitializeTalonDatabase
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_initialize_database \
                --f=$(inputs.gtfFile.path) \
                --g=$(inputs.genomeBuild) \
                --a=$(inputs.annotationVersion) \
                --l=$(inputs.minimumLength) \
                --idprefix=$(inputs.novelPrefix) \
                --5p=$(inputs.cutOff5p) \
                --3p=$(inputs.cutOff3p) \
                --o=$(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: gtfFile
        doc: Gtf annotation containing genes, transcripts, and edges.
        type: File
      - id: genomeBuild
        doc: Name of genome build that the gtf file is based on (ie hg38).
        type: string
      - id: annotationVersion
        doc: Name of supplied annotation (will be used to label data).
        type: string
      - id: minimumLength
        doc: Minimum required transcript length.
        default: 300
        type: int
      - id: novelPrefix
        doc: Prefix for naming novel discoveries in eventual talon runs.
        default: TALON
        type: string
      - id: cutOff5p
        doc: Maximum allowable distance (bp) at the 5' end during annotation.
        default: 500
        type: int
      - id: cutOff3p
        doc: Maximum allowable distance (bp) at the 3' end during annotation.
        default: 300
        type: int
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: memory
        doc: The amount of memory available to the job.
        default: 10G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 60
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: databaseFile
        doc: Talon database.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".db")
  - cwlVersion: v1.2
    id: LabelReads
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_label_reads \
                --f=$(inputs.inputSam.path) \
                --g=$(inputs.referenceGenome.path) \
                --t=$(inputs.threads) \
                --ar=$(inputs.fracaRangeSize) \
                --tmpDir=$(inputs.tmpDir) \
                $(inputs.deleteTmp ? "--deleteTmp" : "") \
                --o=$(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: inputSam
        doc: Sam file of transcripts.
        type: File
      - id: referenceGenome
        doc: Reference genome fasta file.
        type: File
      - id: fracaRangeSize
        doc: Size of post-transcript interval to compute fraction.
        default: 20
        type: int
      - id: tmpDir
        doc: Path to directory for tmp files.
        default: ./tmp_label_reads
        type: string
      - id: deleteTmp
        doc: If set, tmp dir will be removed.
        default: true
        type: boolean
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: threads
        doc: The number of threads to be used.
        default: 4
        type: int
      - id: memory
        doc: The amount of memory available to the job.
        default: 25G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 2880
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: labeledSam
        doc: Sam file with labeled transcripts.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_labeled.sam")
      - id: readLabels
        doc: Tabular file with fraction description per read.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_read_labels.tsv")
  - cwlVersion: v1.2
    id: ReformatGtf
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                talon_reformat_gtf \
                -gtf $(inputs.gtfFile.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: gtfFile
        doc: Gtf annotation containing genes, transcripts, and edges.
        type: File
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: reformattedGtf
        doc: Reformatted gtf file.
        type: File
        outputBinding:
            glob: $(inputs.gtfFile.path)
  - cwlVersion: v1.2
    id: SummarizeDatasets
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                talon_summarize \
                --db $(inputs.databaseFile.path) \
                $(inputs.setVerbose ? "--verbose" : "") \
                --o $(inputs.outputPrefix) \
                $(inputs.datasetGroupsCsv === null ? "" : "--groups " + inputs.datasetGroupsCsv.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: databaseFile
        doc: Talon database.
        type: File
      - id: setVerbose
        doc: Print out the counts in terminal.
        default: false
        type: boolean
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: datasetGroupsCsv
        doc: File of comma-delimited dataset groups to process together.
        type:
          - File
          - 'null'
      - id: memory
        doc: The amount of memory available to the job.
        default: 4G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 50
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: summaryFile
        doc: Tab-delimited file of gene and transcript counts for each dataset.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_talon_summary.tsv")
  - cwlVersion: v1.2
    id: Talon
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                mkdir -p $PWD/tmp #Standard /tmp fills up which makes the SQLite process crash.
                ln -s $PWD/tmp /tmp/sqltmp #Multiprocessing will crash if the absolute path is too long.
                export TMPDIR=/tmp/sqltmp
                printf "" > $(inputs.outputPrefix)/talonConfigFile.csv #File needs to be emptied when task is rerun.
                for file in $(inputs.samFiles.map(function(el) {return el.path}).join(" "))
                do
                    configFileLine="\$(basename ${file%.*}),$(inputs.organism),$(inputs.sequencingPlatform),${file}"
                    echo ${configFileLine} >> $(inputs.outputPrefix)/talonConfigFile.csv
                done
                talon \
                $("--f " + inputs.outputPrefix + "/talonConfigFile.csv") \
                --db $(inputs.databaseFile.path) \
                --build $(inputs.genomeBuild) \
                --threads $(inputs.threads) \
                --cov $(inputs.minimumCoverage) \
                --identity $(inputs.minimumIdentity) \
                $("--o " + inputs.outputPrefix + "/run")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/talon:v5.0_cv1
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
      - id: samFiles
        doc: Input sam files.
        type:
            items: File
            type: array
      - id: organism
        doc: The name of the organism from which the samples originated.
        type: string
      - id: sequencingPlatform
        doc: The sequencing platform used to generate long reads.
        default: PacBio-RS-II
        type: string
      - id: databaseFile
        doc: Talon database. Created using initialize_talon_database.py.
        type: File
      - id: genomeBuild
        doc: Genome build (i.e. hg38) to use.
        type: string
      - id: minimumCoverage
        doc: Minimum alignment coverage in order to use a sam entry.
        default: 0.9
        type: float
      - id: minimumIdentity
        doc: Minimum alignment identity in order to use a sam entry.
        default: 0.8
        type: float
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: threads
        doc: The number of threads to be used.
        default: 4
        type: int
      - id: memory
        doc: The amount of memory available to the job.
        default: 25G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 2880
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/talon:v5.0_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: updatedDatabase
        doc: Updated talon database.
        type: File
        outputBinding:
            glob: $(inputs.databaseFile.path)
      - id: talonLog
        doc: Log file from talon run.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "/run_QC.log")
      - id: talonAnnotation
        doc: Read annotation file from talon run.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "/run_talon_read_annot.tsv")
      - id: talonConfigFile
        doc: The talon configuration file.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "/talonConfigFile.csv")
