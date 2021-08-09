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
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
