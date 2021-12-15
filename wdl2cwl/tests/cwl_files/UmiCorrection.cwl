class: CommandLineTool
id: CorrectUMItools
inputs:
  - id: bam_input
    type: File
  - id: bam_index
    type: File
  - id: docker
    default: quay.io/humancellatlas/secondary-analysis-umitools:0.0.1
    type: string
  - id: output_bam_filename
    default: output.bam
    type: string
  - id: groupout_filename
    default: groupout.tsv
    type: string
  - id: machine_mem_mb
    default: 16000
    type: int
  - id: cpu
    default: 1
    type: int
  - id: preemptible
    default: 3
    type: int
outputs:
  - id: bam_output
    type: File
    outputBinding:
        glob: $(inputs.output_bam_filename)
  - id: group_output
    type: File
    outputBinding:
        glob: $(inputs.groupout_filename)
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/humancellatlas/secondary-analysis-umitools:0.0.1
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
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
  - class: ResourceRequirement
    ramMin: |-
        ${
        var unit = "MiB";
        var value = parseInt(inputs.machine_mem_mb.match(/[0-9]+/g));
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
  - class: ResourceRequirement
    outdirMin: 1024
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
