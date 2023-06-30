cwlVersion: v1.2
id: Mapping
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            pbmm2 align \
            --preset $(inputs.presetOption) \
            $(inputs.sort ? "--sort" : "") \
            -j $(inputs.cores) \
            $(inputs.referenceMMI.path) \
            $(inputs.queryFile.path) \
            --sample $(inputs.sample) \
            $(inputs.sample).align.bam
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/pbmm2:1.3.0--h56fc30b_1
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
  - class: ToolTimeLimit
    timelimit: $(1 + Math.ceil((function(size_of=0){inputs.queryFile.path.forEach(function(element){
        if (element) {size_of += element.size}})}) / 1000^3 * 2000 / inputs.cores)  *
        60)
inputs:
  - id: presetOption
    doc: This option applies multiple options at the same time.
    type: string
  - id: sort
    doc: Sort the output bam file.
    default: true
    type: boolean
  - id: sample
    doc: Name of the sample.
    type: string
  - id: referenceMMI
    doc: MMI file for the reference.
    type: File
  - id: queryFile
    doc: BAM file with reads to align against the reference.
    type: File
  - id: cores
    doc: The number of cores to be used.
    default: 4
    type: int
  - id: memory
    doc: The amount of memory available to the job.
    default: 30G
    type: string
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    type:
      - int
      - 'null'
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/pbmm2:1.3.0--h56fc30b_1
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: outputAlignmentFile
    doc: Mapped bam file.
    type: File
    outputBinding:
        glob: $(inputs.sample + ".align.bam")
  - id: outputIndexFile
    doc: Bam index file.
    type: File
    outputBinding:
        glob: $(inputs.sample + ".align.bam.bai")
