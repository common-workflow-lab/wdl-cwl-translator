class: CommandLineTool
id: Refine
inputs:
  - id: minPolyALength
    default: 20
    type: int
  - id: requirePolyA
    default: false
    type: boolean
  - id: logLevel
    default: WARN
    type: string
  - id: inputBamFile
    type: File
  - id: inputBamIndex
    type: File
  - id: primerFile
    type: File
  - id: outputDir
    type: string
  - id: outputNamePrefix
    type: string
  - id: threads
    default: 2
    type: int
  - id: memory
    default: 2G
    type: string
  - id: timeMinutes
    default: 30
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/isoseq3:3.4.0--0
    type: string
outputs:
  - id: refineBam
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".bam")
  - id: refineBamIndex
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".bam.pbi")
  - id: refineConsensusReadset
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".consensusreadset.xml")
  - id: refineFilterSummary
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".filter_summary.json")
  - id: refineReport
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".report.csv")
  - id: refineStderr
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".stderr.log")
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/isoseq3:3.4.0--0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            set -e
            mkdir -p "$(inputs.outputDir)"
            isoseq3 refine \
            --min-polya-length $(inputs.minPolyALength) \
            $(inputs.requirePolyA ? "--require-polya" : "") \
            --log-level $(inputs.logLevel) \
            --num-threads $(inputs.threads) \
            --log-file "$(inputs.outputDir)/$(inputs.outputNamePrefix).stderr.log" \
            $(inputs.inputBamFile.path) \
            $(inputs.primerFile.path) \
            "$(inputs.outputDir)/$(inputs.outputNamePrefix).bam"
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
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
        return parseInt(memory);
        }
    outdirMin: 1024
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes * 60)
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
