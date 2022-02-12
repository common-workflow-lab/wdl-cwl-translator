class: CommandLineTool
id: Refine
inputs:
  - id: minPolyALength
    doc: Minimum poly(A) tail length.
    default: 20
    type: int
  - id: requirePolyA
    doc: Require fl reads to have a poly(A) tail and remove it.
    default: false
    type: boolean
  - id: logLevel
    doc: 'Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).'
    default: WARN
    type: string
  - id: inputBamFile
    doc: Bam input file.
    type: File
  - id: inputBamIndex
    doc: Index for the Bam input file.
    type: File
  - id: primerFile
    doc: Barcode/primer fasta file.
    type: File
  - id: outputDir
    doc: Output directory path.
    type: string
  - id: outputNamePrefix
    doc: Basename of the output files.
    type: string
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
    default: 30
    type: int
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/isoseq3:3.4.0--0
    type: string
outputs:
  - id: refineBam
    doc: Filtered reads output file.
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".bam")
  - id: refineBamIndex
    doc: Index of filtered reads output file.
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".bam.pbi")
  - id: refineConsensusReadset
    doc: Refine consensus readset xml file.
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".consensusreadset.xml")
  - id: refineFilterSummary
    doc: Refine summary file.
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".filter_summary.json")
  - id: refineReport
    doc: Refine report file.
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".report.csv")
  - id: refineStderr
    doc: Refine stderr log file.
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.outputNamePrefix + ".stderr.log")
requirements:
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
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/isoseq3:3.4.0--0
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
