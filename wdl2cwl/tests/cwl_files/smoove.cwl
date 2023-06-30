cwlVersion: v1.2
id: Call
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            set -e
            mkdir -p $(inputs.outputDir)
            smoove call \
            --outdir $(inputs.outputDir) \
            --name $(inputs.sample) \
            --fasta $(inputs.referenceFasta.path) \
            --removepr \
            --genotype \
            $(inputs.bamFile.path)
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/smoove:0.2.5--0
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
  - id: bamFile
    doc: The bam file to process.
    type: File
  - id: bamIndex
    doc: The index of the bam file.
    type: File
  - id: referenceFasta
    doc: The reference fasta file also used for mapping.
    type: File
  - id: referenceFastaFai
    doc: Fasta index (.fai) file of the reference.
    type: File
  - id: sample
    doc: The name of the sample.
    type: string
  - id: outputDir
    doc: The location the output VCF file should be written.
    default: ./smoove
    type: string
  - id: memory
    doc: The memory required to run the programs.
    default: 15G
    type: string
  - id: timeMinutes
    doc: The maximum duration (in minutes) the tool is allowed to run.
    default: 1440
    type: int
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/smoove:0.2.5--0
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: smooveVcf
    doc: Calls of structural variants in VCF file.
    type: File
    outputBinding:
        glob: $(inputs.outputDir + "/" + inputs.sample + "-smoove.genotyped.vcf.gz")
