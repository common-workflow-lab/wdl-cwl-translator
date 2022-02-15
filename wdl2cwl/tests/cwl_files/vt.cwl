cwlVersion: v1.2
id: Normalize
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            set -eo pipefail
            mkdir -p "\$(dirname $(inputs.outputPath))"
            vt normalize $(inputs.inputVCF.path) \
            -r $(inputs.referenceFasta.path) \
            $(inputs.ignoreMaskedRef ? "-m " : "") \
            | vt decompose -s - -o $(inputs.outputPath)
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/vt:0.57721--hdf88d34_2
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
        else throw "Unknown units: " + unit;
        return parseInt(memory);
        }
    outdirMin: 1024
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes * 60)
inputs:
  - id: inputVCF
    doc: The VCF file to process.
    type: File
  - id: inputVCFIndex
    doc: The index of the VCF file to be processed.
    type: File
  - id: referenceFasta
    doc: The reference fasta file which was also used for mapping.
    type: File
  - id: referenceFastaFai
    doc: The index for the reference fasta file.
    type: File
  - id: ignoreMaskedRef
    doc: Warns but does not exit when REF is inconsistent with masked reference sequence
        for non SNPs.
    default: false
    type: boolean
  - id: outputPath
    doc: The location the output VCF file should be written.
    default: ./vt/normalized_decomposed.vcf
    type: string
  - id: memory
    doc: The memory required to run the programs.
    default: 4G
    type: string
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    default: 30
    type: int
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/vt:0.57721--hdf88d34_2
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: outputVcf
    doc: Normalized & decomposed VCF file.
    type: File
    outputBinding:
        glob: $(inputs.outputPath)
