class: CommandLineTool
id: TrimAdapters
inputs:
  - id: fastq1_input_files
    type:
        items: File
        type: array
  - id: fastq2_input_files
    type:
        items: File
        type: array
  - id: adapter_list
    type: File
  - id: input_ids
    type:
        items: string
        type: array
  - id: docker
    default: quay.io/humancellatlas/snss2-trim-adapters:0.1.0
    type: string
  - id: machine_mem_mb
    default: 8250
    type: int
  - id: cpu
    default: 1
    type: int
  - id: disk
    type:
      - int
      - 'null'
  - id: preemptible
    default: 3
    type: int
outputs:
  - id: trimmed_fastq1_files
    type:
        items: File
        type: array
    outputBinding:
        glob: $("*trimmed_R1.fastq.gz")
  - id: trimmed_fastq2_files
    type:
        items: File
        type: array
    outputBinding:
        glob: $("*trimmed_R2.fastq.gz")
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/humancellatlas/snss2-trim-adapters:0.1.0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e

            declare -a fastq1_files=($(inputs.fastq1_input_files.map(function(el) {return el.path}).join(" ")))
            declare -a fastq2_files=($(inputs.fastq2_input_files.map(function(el) {return el.path}).join(" ")))
            declare -a output_prefix=($(inputs.input_ids.join(" ")))
            for ((i=0; i<${#fastq1_files[@]}; ++i));
              do
                fastq1=${fastq1_files[$i]}
                fastq2=${fastq2_files[$i]}

                fastq-mcf \
                   -C 200000 $(inputs.adapter_list.path) \
                   $fastq1 \
                   $fastq2 \
                   -o "${output_prefix[$i]}.trimmed_R1.fastq.gz" \
                   -o "${output_prefix[$i]}.trimmed_R2.fastq.gz"
              done;
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
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
        return parseInt(memory);
        }
    outdirMin: $('local-disk ' + $(Math.ceil(2*(inputs.fastq1_input_files.size / 1024^3
        + inputs.fastq2_input_files.size / 1024^3)) + 10) + ' HDD')
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
