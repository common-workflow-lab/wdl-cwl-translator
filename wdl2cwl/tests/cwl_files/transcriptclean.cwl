cwlVersion: v1.2
$graph:
  - cwlVersion: v1.2
    id: GetSJsFromGtf
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                get_SJs_from_gtf \
                --f=$(inputs.gtfFile.path) \
                --g=$(inputs.genomeFile.path) \
                --minIntronSize=$(inputs.minIntronSize) \
                $("--o=" + inputs.outputPrefix + ".tsv")
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/transcriptclean:v2.0.2_cv1
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
        doc: Input gtf file.
        type: File
      - id: genomeFile
        doc: Reference genome.
        type: File
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: minIntronSize
        doc: Minimum size of intron to consider a junction.
        default: 21
        type: int
      - id: memory
        doc: The amount of memory available to the job.
        default: 8G
        type: string
      - id: timeMinutes
        doc: The maximum amount of time the job will run in minutes.
        default: 30
        type: int
      - id: dockerImage
        doc: The docker image used for this task. Changing this may result in errors
            which the developers may choose not to address.
        default: biocontainers/transcriptclean:v2.0.2_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: spliceJunctionFile
        doc: Extracted splice junctions.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".tsv")
  - cwlVersion: v1.2
    id: GetTranscriptCleanStats
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                get_TranscriptClean_stats \
                $(inputs.inputSam.path) \
                $(inputs.outputPrefix)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/transcriptclean:v2.0.2_cv1
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
      - id: inputSam
        doc: Output sam file from transcriptclean.
        type: File
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
        default: biocontainers/transcriptclean:v2.0.2_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: statsFile
        doc: Summary stats from transcriptclean run.
        type: stdout
  - cwlVersion: v1.2
    id: TranscriptClean
    class: CommandLineTool
    requirements:
      - class: InitialWorkDirRequirement
        listing:
          - entryname: script.bash
            entry: |4

                set -e
                mkdir -p "\$(dirname $(inputs.outputPrefix))"
                TranscriptClean \
                -s $(inputs.inputSam.path) \
                -g $(inputs.referenceGenome.path) \
                -t $(inputs.cores) \
                --maxLenIndel=$(inputs.maxLenIndel) \
                --maxSJOffset=$(inputs.maxSJOffset) \
                -o $(inputs.outputPrefix) \
                $(inputs.correctMismatches ? "-m true" : "-m false") \
                $(inputs.correctIndels ? "-i true" : "-i false") \
                $(inputs.correctSJs ? "--correctSJs=true" : "--correctSJs=false") \
                $(inputs.dryRun ? "--dryRun" : "") \
                $(inputs.primaryOnly ? "--primaryOnly" : "") \
                $(inputs.canonOnly ? "--canonOnly" : "") \
                --bufferSize=$(inputs.bufferSize) \
                $(inputs.deleteTmp ? "--deleteTmp" : "") \
                $(inputs.spliceJunctionAnnotation === null ? "" : "-j " + inputs.spliceJunctionAnnotation.path) \
                $(inputs.variantFile === null ? "" : "-v " + inputs.variantFile.path)
      - class: InlineJavascriptRequirement
      - class: NetworkAccess
        networkAccess: true
    hints:
      - class: DockerRequirement
        dockerPull: biocontainers/transcriptclean:v2.0.2_cv1
      - class: ResourceRequirement
        coresMin: $(inputs.cores)
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
    inputs:
      - id: inputSam
        doc: Input sam file containing transcripts to correct.
        type: File
      - id: referenceGenome
        doc: Reference genome fasta file.
        type: File
      - id: maxLenIndel
        doc: Maximum size indel to correct.
        default: 5
        type: int
      - id: maxSJOffset
        doc: Maximum distance from annotated splice junction to correct.
        default: 5
        type: int
      - id: outputPrefix
        doc: Output directory path + output file prefix.
        type: string
      - id: correctMismatches
        doc: Set this to make transcriptclean correct mismatches.
        default: true
        type: boolean
      - id: correctIndels
        doc: Set this to make transcriptclean correct indels.
        default: true
        type: boolean
      - id: correctSJs
        doc: Set this to make transcriptclean correct splice junctions.
        default: true
        type: boolean
      - id: dryRun
        doc: Transcriptclean will read in the data but don't do any correction.
        default: false
        type: boolean
      - id: primaryOnly
        doc: Only output primary mappings of transcripts.
        default: false
        type: boolean
      - id: canonOnly
        doc: Only output canonical transcripts and transcript containing annotated
            noncanonical junctions.
        default: false
        type: boolean
      - id: bufferSize
        doc: Number of lines to output to file at once by each thread during run.
        default: 100
        type: int
      - id: deleteTmp
        doc: The temporary directory generated by transcriptclean will be removed.
        default: true
        type: boolean
      - id: spliceJunctionAnnotation
        doc: Splice junction file.
        type:
          - File
          - 'null'
      - id: variantFile
        doc: Vcf formatted file of variants.
        type:
          - File
          - 'null'
      - id: cores
        doc: The number of cores to be used.
        default: 1
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
        default: biocontainers/transcriptclean:v2.0.2_cv1
        type: string
    baseCommand:
      - bash
      - script.bash
    outputs:
      - id: fastaFile
        doc: Fasta file containing corrected reads.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.fa")
      - id: logFile
        doc: Log file of transcriptclean run.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.log")
      - id: outputSam
        doc: Sam file containing corrected aligned reads.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.sam")
      - id: logFileTE
        doc: TE log file of transcriptclean run.
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.TE.log")
