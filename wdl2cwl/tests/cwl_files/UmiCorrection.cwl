cwlVersion: v1.2
id: CorrectUMItools
class: CommandLineTool
doc: Marks duplicates using umitools group specifically for single-cell experiments
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4+

             set -e

             mv $(inputs.bam_input.path) input.bam
             mv $(inputs.bam_index.path) input.bam.bai

             touch input.bam
             touch input.bam.bai

             umi_tools group \
                 -I input.bam \
                 -L outlog.txt \
                 -E outerr.txt \
                 -S duplicate_marked.bam \
                 --output-bam \
                 --extract-umi-method=tag \
                 --umi-tag UR \
                 --method directional \
                 --per-gene \
                 --per-cell \
                 --cell-tag CB \
                 --gene-tag GE \
                 --no-sort-output \
                 --group-out $(inputs.groupout_filename) \
                 --umi-group-tag UB

            getUntaggedReads --in-bam-file input.bam --out-bam-file untagged.bam

            rm input.bam input.bam.bai
            samtools cat -o $(inputs.output_bam_filename) duplicate_marked.bam untagged.bam

  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/humancellatlas/secondary-analysis-umitools:0.0.1
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: |-
        ${
        var unit = "MiB";
        var value = parseInt(`${inputs.machine_mem_mb}`.match(/[0-9]+/g));
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
    outdirMin: $((Math.ceil((function(size_of=0){inputs.bam_input.forEach(function(element){
        if (element) {size_of += element.size}})}) / 1024^3 * 6)  + 50) * 1024)
inputs:
  - id: bam_input
    doc: Aligned and sorted bam
    type: File
  - id: bam_index
    type: File
  - id: docker
    doc: (optional) the docker image containing the runtime environment for this task
    default: quay.io/humancellatlas/secondary-analysis-umitools:0.0.1
    type: string
  - id: output_bam_filename
    default: output.bam
    type: string
  - id: groupout_filename
    default: groupout.tsv
    type: string
  - id: machine_mem_mb
    doc: (optional) the amount of memory (MiB) to provision for this task
    default: 16000
    type: int
  - id: cpu
    doc: (optional) the number of cpus to provision for this task
    default: 1
    type: int
  - id: disk
    doc: (optional) the amount of disk space (GiB) to provision for this task
    type:
      - int
      - 'null'
  - id: preemptible
    doc: (optional) if non-zero, request a pre-emptible instance and allow for this
        number of preemptions before running the task on a non preemptible machine
    default: 3
    type: int
baseCommand:
  - bash
  - script.bash
outputs:
  - id: bam_output
    type: File
    outputBinding:
        glob: $(inputs.output_bam_filename)
  - id: group_output
    type: File
    outputBinding:
        glob: $(inputs.groupout_filename)
