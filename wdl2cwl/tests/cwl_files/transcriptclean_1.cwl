class: CommandLineTool
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
  - class: DockerRequirement
    dockerPull: biocontainers/transcriptclean:v2.0.2_cv1
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p "\$(dirname $(inputs.outputPrefix))"
            get_TranscriptClean_stats \
            $(inputs.inputSam.path) \
            $(inputs.outputPrefix)
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: |-
        ${
        var unit = inputs["memory"].match(/[a-zA-Z]+/g).join("");
        var value = parseInt(inputs["memory"].match(/[0-9]+/g));
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
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes* 60)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
