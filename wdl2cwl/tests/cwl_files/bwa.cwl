cwlVersion: v1.2
id: Mem
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            set -e
            mkdir -p "\$(dirname $(inputs.outputPrefix))"
            bwa mem \
              -t $(inputs.threads) \
              $(inputs.useSoftclippingForSupplementary ? "-Y" : "") \
              $(inputs.readgroup === null ? "" : "-R '" + inputs.readgroup)$(inputs.readgroup === null ? "" : "'") \
              $(inputs.bwaIndex.fastaFile.path) \
              $(inputs.read1.path) \
              $(inputs.read2 === null ? "" : inputs.read2.path) \
              2> $(inputs.outputPrefix).log.bwamem | \
              $(inputs.usePostalt ? "" : "#") bwa-postalt.js -p $(inputs.outputPrefix).hla $(inputs.bwaIndex.fastaFile.path)$(inputs.sixtyFour ? ".64.alt" : ".alt") | \
              samtools sort \
              $("-@ " + inputs.totalSortThreads) \
              -m $(inputs.sortMemoryPerThreadGb)G \
              -l $(inputs.compressionLevel) \
              - \
              -o $(inputs.outputPrefix).aln.bam
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: |-
        ${
        var unit = "G";
        var value = parseInt(`${[inputs.memoryGb, 10 + Math.ceil((function(size_of=0){inputs.bwaIndex.indexFiles.forEach(function(element){ if (element) {size_of += element.size}})}) / 1000^3 * 2)  + inputs.sortMemoryPerThreadGb * [inputs.sortThreads, inputs.threads === 1 ? 1 : 1 + Math.ceil(inputs.threads / 4.0) ].find(function(element) { return element !== null }) ].find(function(element) { return element !== null }) }`.match(/[0-9]+/g));
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
    timelimit: '$(10 + Math.ceil((function(size_of=0){[inputs.read1.path, inputs.read2
        === null ? "" : inputs.read2.path].forEach(function(element){ if (element)
        {size_of += element.size}})}) / 1000^3 * 300 / inputs.threads)  * 60)'
inputs:
  - id: read1
    doc: The first-end fastq file.
    type: File
  - id: read2
    doc: The second-end fastq file.
    type:
      - File
      - 'null'
  - id: bwaIndex
    doc: The BWA index, including (optionally) a .alt file.
    type:
        name: BwaIndex
        fields:
          - name: fastaFile
            type: File
          - name: indexFiles
            type:
                items: File
                type: array
        type: record
  - id: outputPrefix
    doc: The prefix of the output files, including any parent directories.
    type: string
  - id: sixtyFour
    doc: Whether or not the index uses the '.64' suffixes.
    default: false
    type: boolean
  - id: usePostalt
    doc: Whether to use the postalt script from bwa kit.
    default: false
    type: boolean
  - id: useSoftclippingForSupplementary
    default: false
    type: boolean
  - id: sortMemoryPerThreadGb
    doc: The amount of memory for each sorting thread in gigabytes.
    default: 2
    type: int
  - id: compressionLevel
    doc: The compression level of the output BAM.
    default: 1
    type: int
  - id: readgroup
    doc: A readgroup identifier.
    type:
      - string
      - 'null'
  - id: sortThreads
    doc: The number of threads to use for sorting.
    type:
      - int
      - 'null'
  - id: threads
    doc: The number of threads to use for alignment.
    default: 4
    type: int
  - id: memoryGb
    doc: The amount of memory this job will use in gigabytes.
    type:
      - int
      - 'null'
  - id: timeMinutes
    doc: The maximum amount of time the job will run in minutes.
    type:
      - int
      - 'null'
  - id: dockerImage
    doc: The docker image used for this task. Changing this may result in errors which
        the developers may choose not to address.
    default: quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0
    type: string
baseCommand:
  - bash
  - script.bash
outputs:
  - id: outputBam
    doc: The produced BAM file.
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix + ".aln.bam")
  - id: outputHla
    doc: The produced HLA file.
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputPrefix + ".hla")
