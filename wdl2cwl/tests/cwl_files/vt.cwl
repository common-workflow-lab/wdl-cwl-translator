class: CommandLineTool
id: Normalize
inputs:
  - id: inputVCF
    type: File
  - id: inputVCFIndex
    type: File
  - id: referenceFasta
    type: File
  - id: referenceFastaFai
    type: File
  - id: ignoreMaskedRef
    default: false
    type: boolean
  - id: outputPath
    default: ./vt/normalized_decomposed.vcf
    type: string
  - id: memory
    default: 4G
    type: string
  - id: timeMinutes
    default: 30
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/vt:0.57721--hdf88d34_2
    type: string
outputs:
  - id: outputVcf
    type: File
    outputBinding:
        glob: $(inputs.outputPath)
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/vt:0.57721--hdf88d34_2
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
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
  - class: ResourceRequirement
    ramMin: |-
        ${
        var unit = inputs.memory.match(/[a-zA-Z]+/g).join("");
        var value = parseInt(inputs.memory.match(/[0-9]+/g));
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
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes * 60)
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
