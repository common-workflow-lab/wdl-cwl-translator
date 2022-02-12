cwlVersion: v1.2
id: VarDict
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            set -e -o pipefail
            export JAVA_OPTS="-Xmx$(inputs.javaXmx) -XX:ParallelGCThreads=1"
            vardict-java \
            -th  $(inputs.threads) \
            -G $(inputs.referenceFasta.path) \
            -N $(inputs.tumorSampleName) \
            -b "$(inputs.tumorBam.path)$(inputs.normalBam === null ? "" : "|" + inputs.normalBam.path)" \
            $(inputs.normalBam === null ? "-z" : "") \
            -c $(inputs.chromosomeColumn) \
            -S $(inputs.startColumn) \
            -E $(inputs.endColumn) \
            -g $(inputs.geneColumn) \
            $(inputs.bedFile.path) | \
            $(inputs.normalBam === null ? "teststrandbias.R" : "testsomatic.R") | \
            $(inputs.normalBam === null ? "var2vcf_valid.pl" : "var2vcf_paired.pl") \
            -N "$(inputs.tumorSampleName)$(inputs.normalSampleName === null ? "" : "|" + inputs.normalSampleName)" \
            $(inputs.normalBam === null ? "-E" : "") \
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
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/vardict-java:1.5.8--1
  - class: ResourceRequirement
    coresMin: $(inputs.threads + 2)
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
        return parseInt(memory);
        }
    outdirMin: 1024
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes * 60)
inputs:
  - id: tumorSampleName
    doc: The name of the tumor/case sample.
    type: string
  - id: tumorBam
    doc: The tumor/case sample's BAM file.
    type: File
  - id: tumorBamIndex
    doc: The index for the tumor/case sample's BAM file.
    type: File
  - id: referenceFasta
    doc: The reference fasta file.
    type: File
  - id: referenceFastaFai
    doc: The index for the reference fasta file.
    type: File
  - id: bedFile
    doc: A bed file describing the regions to operate on. These regions must be below
        1e6 bases in size.
    type: File
  - id: outputVcf
    doc: The location to write the output VCF file to.
    type: string
  - id: outputCandidateSomaticOnly
    doc: Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-M` flag.
    default: true
    type: boolean
  - id: outputAllVariantsAtSamePosition
    doc: Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-A` flag.
    default: true
    type: boolean
  - id: mappingQuality
    doc: Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-Q` option.
    default: 20.0
    type: float
  - id: minimumTotalDepth
    doc: Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-d` option.
    default: 8
    type: int
  - id: minimumVariantDepth
    doc: Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-v` option.
    default: 4
    type: int
  - id: minimumAlleleFrequency
    doc: Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-f` option.
    default: 0.02
    type: float
  - id: chromosomeColumn
    doc: Equivalent to vardict-java's `-c` option.
    default: 1
    type: int
  - id: startColumn
    doc: Equivalent to vardict-java's `-S` option.
    default: 2
    type: int
  - id: endColumn
    doc: Equivalent to vardict-java's `-E` option.
    default: 3
    type: int
  - id: geneColumn
    doc: Equivalent to vardict-java's `-g` option.
    default: 4
    type: int
  - id: normalSampleName
    doc: The name of the normal/control sample.
    type:
      - string
      - 'null'
  - id: normalBam
    doc: The normal/control sample's BAM file.
    type:
      - File
      - 'null'
  - id: normalBamIndex
    doc: The normal/control sample's BAM file.
    type:
      - File
      - 'null'
  - id: javaXmx
    doc: The maximum memory available to the program. Should be lower than `memory`
        to accommodate JVM overhead.
    default: 16G
    type: string
  - id: threads
    doc: The number of threads to use.
    default: 1
    type: int
  - id: memory
    doc: The amount of memory this job will use.
    default: 18G
    type: string
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    default: 300
    type: int
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/vardict-java:1.5.8--1
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: vcfFile
    doc: Output VCF file.
    type: File
    outputBinding:
        glob: $(inputs.outputVcf)
