class: CommandLineTool
id: Mem
inputs:
  - id: read1
    type: File
  - id: read2
    type:
      - File
      - 'null'
  - id: outputPrefix
    type: string
  - id: readgroup
    type:
      - string
      - 'null'
  - id: sortThreads
    type:
      - int
      - 'null'
  - id: memoryGb
    type:
      - int
      - 'null'
  - id: estimatedMemoryGb
    type: int
  - id: sixtyFour
    default: false
    type: boolean
  - id: usePostalt
    default: false
    type: boolean
  - id: sortMemoryPerThreadGb
    default: 2
    type: int
  - id: compressionLevel
    default: 1
    type: int
  - id: threads
    default: 4
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0
    type: string
outputs:
  - id: outputBam
    type: File
    outputBinding:
        glob: $(inputs.outputPrefix).aln.bam
  - id: outputHla
    type:
      - File
      - 'null'
    outputBinding:
        glob: $(inputs.outputPrefix).hla
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4

            set -e
            mkdir -p "\$(dirname $(inputs.outputPrefix))"
            bwa mem \
              -t $(inputs.threads) \
              $(inputs.readgroup === null ? "" : "-R '" + inputs.readgroup )$(inputs.readgroup === null  ? "" : "'") \
              $(inputs.bwaIndex.fastaFile) \
              $(inputs.read1.path) \
              $(inputs.read2 === null ? "" : inputs.read2.path) \
              2> $(inputs.outputPrefix).log.bwamem | \
              $(inputs.usePostalt ? "" : "#") bwa-postalt.js -p $(inputs.outputPrefix).hla $(inputs.bwaIndex.fastaFile)$(inputs.sixtyFour ? ".64.alt" : ".alt") | \
              samtools sort \
              -m $(inputs.sortMemoryPerThreadGb)G \
              -l $(inputs.compressionLevel) \
              - \
              -o $(inputs.outputPrefix).aln.bam
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    ramMin: |-
        ${
        var unit = "G;
        var value = (function() {for (const elem of [inputs.memoryGb,inputs.estimatedMemoryGb]) if (elem != null) return elem}) ();
        if (value == undefined) throw "error! array contains only null values or it's empty";
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
    coresMin: $(inputs.threads)
cwlVersion: v1.2
baseCommand:
  - bash
  - example.sh
