cwlVersion: v1.2
$graph:
  - class: CommandLineTool
    id: GetSJsFromGtf
    inputs:
      - id: gtfFile
        type: File
      - id: genomeFile
        type: File
      - id: outputPrefix
        type: string
      - id: minIntronSize
        default: 21
        type: int
      - id: memory
        default: 8G
        type: string
      - id: timeMinutes
        default: 30
        type: int
      - id: dockerImage
        default: biocontainers/transcriptclean:v2.0.2_cv1
        type: string
    outputs:
      - id: spliceJunctionFile
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + ".tsv")
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
                --o= $(inputs.outputPrefix) + ".tsv"
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: GetTranscriptCleanStats
    inputs:
      - id: inputSam
        type: File
      - id: outputPrefix
        type: string
      - id: memory
        default: 4G
        type: string
      - id: timeMinutes
        default: 30
        type: int
      - id: dockerImage
        default: biocontainers/transcriptclean:v2.0.2_cv1
        type: string
    outputs:
      - id: statsFile
        type: stdout
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
            return parseInt(memory);
            }
        outdirMin: 1024
      - class: ToolTimeLimit
        timelimit: $(inputs.timeMinutes * 60)
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
  - class: CommandLineTool
    id: TranscriptClean
    inputs:
      - id: inputSam
        type: File
      - id: referenceGenome
        type: File
      - id: maxLenIndel
        default: 5
        type: int
      - id: maxSJOffset
        default: 5
        type: int
      - id: outputPrefix
        type: string
      - id: correctMismatches
        default: true
        type: boolean
      - id: correctIndels
        default: true
        type: boolean
      - id: correctSJs
        default: true
        type: boolean
      - id: dryRun
        default: false
        type: boolean
      - id: primaryOnly
        default: false
        type: boolean
      - id: canonOnly
        default: false
        type: boolean
      - id: bufferSize
        default: 100
        type: int
      - id: deleteTmp
        default: true
        type: boolean
      - id: spliceJunctionAnnotation
        type:
          - File
          - 'null'
      - id: variantFile
        type:
          - File
          - 'null'
      - id: cores
        default: 1
        type: int
      - id: memory
        default: 25G
        type: string
      - id: timeMinutes
        default: 2880
        type: int
      - id: dockerImage
        default: biocontainers/transcriptclean:v2.0.2_cv1
        type: string
    outputs:
      - id: fastaFile
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.fa")
      - id: logFile
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.log")
      - id: outputSam
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.sam")
      - id: logFileTE
        type: File
        outputBinding:
            glob: $(inputs.outputPrefix + "_clean.TE.log")
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
            return parseInt(memory);
            }
        outdirMin: 1024
    cwlVersion: v1.2
    baseCommand:
      - bash
      - script.bash
