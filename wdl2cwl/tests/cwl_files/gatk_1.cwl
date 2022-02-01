class: CommandLineTool
id: HaplotypeCaller
inputs:
  - id: inputBams
    type:
        items: File
        type: array
  - id: inputBamsIndex
    type:
        items: File
        type: array
  - id: outputPath
    type: string
  - id: referenceFasta
    type: File
  - id: referenceFastaIndex
    type: File
  - id: referenceFastaDict
    type: File
  - id: gvcf
    default: false
    type: boolean
  - id: dontUseSoftClippedBases
    default: false
    type: boolean
  - id: intervalList
    type:
      - items: File
        type: array
      - 'null'
  - id: excludeIntervalList
    type:
      - items: File
        type: array
      - 'null'
  - id: contamination
    type:
      - float
      - 'null'
  - id: dbsnpVCF
    type:
      - File
      - 'null'
  - id: dbsnpVCFIndex
    type:
      - File
      - 'null'
  - id: pedigree
    type:
      - File
      - 'null'
  - id: ploidy
    type:
      - int
      - 'null'
  - id: outputMode
    type:
      - string
      - 'null'
  - id: standardMinConfidenceThresholdForCalling
    type:
      - float
      - 'null'
  - id: javaXmxMb
    default: 4096
    type: int
  - id: timeMinutes
    default: 400
    type: int
  - id: dockerImage
    default: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
    type: string
outputs:
  - id: outputVCF
    type: File
    outputBinding:
        glob: $(inputs.outputPath)
  - id: outputVCFIndex
    type: File
    outputBinding:
        glob: $(inputs.outputPath + ".tbi")
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4+

            set -e
            mkdir -p "\$(dirname $(inputs.outputPath))"
            mkdir wd
            for FILE in $(inputs.inputBams.map(function(el) {return el.path}).join(" ")); do ln -s $FILE wd/\$(inputBams $FILE) ; done
            for FILE in $(inputs.inputBamsIndex.map(function(el) {return el.path}).join(" ")); do ln -s $FILE wd/\$(inputBamsIndex $FILE) ; done
            mkdir wd2
            ln -s $(inputs.referenceFasta.path) wd2/\$(basename $(inputs.referenceFasta.path))
            ln -s $(inputs.referenceFastaDict.path) wd2/\$(basename $(inputs.referenceFastaDict.path))
            ln -s $(inputs.referenceFastaIndex.path) wd2/\$(basename $(inputs.referenceFastaIndex.path))
            gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
            HaplotypeCaller \
            -R wd2/\$(basename $(inputs.referenceFasta.path)) \
            -O $(inputs.outputPath) \
            (for FILE in $(inputs.inputBams.map(function(el) {return el.path}).join(" ")); do echo -- "-I wd/"\$(basename $FILE); done)
            $(inputs.intervalList === null ? "" : "-L") $(inputs.intervalList.map(function(el) {return el.path}).join(" -L ")) \
            $(inputs.excludeIntervalList === null ? "" : "-XL") $(inputs.excludeIntervalList.map(function(el) {return el.path}).join(" -XL ")) \
            $(inputs.dontUseSoftClippedBases ? "--dont-use-soft-clipped-bases" : "") \

  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    outdirMin: 1024
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes * 60)
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
