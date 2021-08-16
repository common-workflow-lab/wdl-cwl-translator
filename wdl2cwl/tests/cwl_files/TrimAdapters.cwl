class: CommandLineTool
id: TrimAdapters
inputs:
  - id: fastq1
    type: File
  - id: fastq2
    default: ''
    type:
      - File
      - 'null'
  - id: adapter_list
    type: File
  - id: docker
    default: quay.io/humancellatlas/snss2-trim-adapters:0.1.0
    type: string
  - id: machine_mem_mb
    default: 8250
    type: int
  - id: cpu
    default: 1
    type: int
  - id: preemptible
    default: 3
    type: int
outputs:
  - id: trimmed_fastq1
    type: File
    outputBinding:
        glob: fastq_R1.trimmed.fastq.gz
  - id: trimmed_fastq2
    type: File
    outputBinding:
        glob: fastq_R2.trimmed.fastq.gz
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/humancellatlas/snss2-trim-adapters:0.1.0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e

            fastq-mcf \
               -C 200000 $(inputs.adapter_list.path) \
               $(inputs.fastq1.path) \
               $(inputs.fastq2) \
               -o fastq_R1.trimmed.fastq.gz \
               -o fastq_R2.trimmed.fastq.gz
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
