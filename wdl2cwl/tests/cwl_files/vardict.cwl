class: CommandLineTool
id: VarDict
inputs:
  - id: tumorSampleName
    type: string
  - id: tumorBam
    type: File
  - id: tumorBamIndex
    type: File
  - id: referenceFasta
    type: File
  - id: referenceFastaFai
    type: File
  - id: bedFile
    type: File
  - id: outputVcf
    type: string
  - id: normalSampleName
    type:
      - string
      - 'null'
  - id: normalBam
    type:
      - File
      - 'null'
  - id: normalBamIndex
    type:
      - File
      - 'null'
  - id: outputCandidateSomaticOnly
    default: true
    type: boolean
  - id: outputAllVariantsAtSamePosition
    default: true
    type: boolean
  - id: mappingQuality
    default: 20.0
    type: float
  - id: minimumTotalDepth
    default: 8
    type: int
  - id: minimumVariantDepth
    default: 4
    type: int
  - id: minimumAlleleFrequency
    default: 0.02
    type: float
  - id: chromosomeColumn
    default: 1
    type: int
  - id: startColumn
    default: 2
    type: int
  - id: endColumn
    default: 3
    type: int
  - id: geneColumn
    default: 4
    type: int
  - id: javaXmx
    default: 16G
    type: string
  - id: threads
    default: 1
    type: int
  - id: memory
    default: 18G
    type: string
  - id: timeMinutes
    default: 300
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/vardict-java:1.5.8--1
    type: string
outputs:
  - id: vcfFile
    type: File
    outputBinding:
        glob: $(inputs.outputVcf)
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/vardict-java:1.5.8--1
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e -o pipefail
            export JAVA_OPTS="-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1"
            vardict-java \
            -th $(inputs.threads) \
            -G $(inputs.referenceFasta.path) \
            -N $(inputs.tumorSampleName) \
            -b "$(inputs.tumorBam.path)$(inputs.normalBam === null ? "" : "|" + inputs.normalBam.path )" \
            $(inputs.normalBam === null  ? "-z" : "") \
            -c $(inputs.chromosomeColumn) \
            -S $(inputs.startColumn) \
            -E $(inputs.endColumn) \
            -g $(inputs.geneColumn) \
            $(inputs.bedFile.path) | \
            $(inputs.normalBam === null  ? "teststrandbias.R" : "testsomatic.R") | \
            $(inputs.normalBam === null  ? "var2vcf_valid.pl" : "var2vcf_paired.pl") \
            -N "$(inputs.tumorSampleName)$(inputs.normalSampleName === null ? "" : "|" + inputs.normalSampleName )" \
            $(inputs.normalBam === null  ? "-E" : "") \
            $(inputs.outputCandidateSomaticOnly ? "-M" : "") \
            $(inputs.outputAllVariantsAtSamePosition ? "-A" : "") \
            -Q $(inputs.mappingQuality) \
            -d $(inputs.minimumTotalDepth) \
            -v $(inputs.minimumVariantDepth) \
            -f $(inputs.minimumAlleleFrequency) \
            > $(inputs.outputVcf)
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
  - class: ResourceRequirement
    coresMin: $(inputs.threads + 2)
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
