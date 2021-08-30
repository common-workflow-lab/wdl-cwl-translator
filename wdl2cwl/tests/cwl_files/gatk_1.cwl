class: CommandLineTool
id: HaplotypeCaller
inputs:
  - id: inputBams
    type:
      - items: File
        type: array
  - id: inputBamsIndex
    type:
      - items: File
        type: array
  - id: outputPath
    type: string
  - id: referenceFasta
    type: File
  - id: referenceFastaIndex
    type: File
  - id: referenceFastaDict
    type: File
  - id: intervalList
    default: []
    type:
      - items:
          - File
          - 'null'
        type: array
  - id: excludeIntervalList
    default: []
    type:
      - items:
          - File
          - 'null'
        type: array
  - id: gvcf
    default: false
    type: boolean
  - id: dontUseSoftClippedBases
    default: false
    type: boolean
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
        glob: $(inputs.outputPath).tbi
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0
  - class: InitialWorkDirRequirement
    listing:
      - entryname: example.sh
        entry: |4+

            set -e
            mkdir -p "\$(dirname $(inputs.outputPath))"
            mkdir wd
            for FILE in $(" ".join(inputs["inputBams"].map(function(el) { return el.path}))); do ln -s $FILE wd/\$(inputBams $FILE) ; done
            for FILE in $(" ".join(inputs["inputBamsIndex"].map(function(el) { return el.path}))); do ln -s $FILE wd/\$(inputBamsIndex $FILE) ; done
            mkdir wd2
            ln -s $(inputs.referenceFasta.path) wd2/\$(basename $(inputs.referenceFasta.path))
            ln -s $(inputs.referenceFastaDict.path) wd2/\$(basename $(inputs.referenceFastaDict.path))
            ln -s $(inputs.referenceFastaIndex.path) wd2/\$(basename $(inputs.referenceFastaIndex.path))
            gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
            HaplotypeCaller \
            -R wd2/\$(basename $(inputs.referenceFasta.path)) \
            -O $(inputs.outputPath) \
            (for FILE in $(" ".join(inputs["inputBams"].map(function(el) { return el.path}))); do echo -- "-I wd/"\$(basename $FILE); done)
            $(inputs["intervalList"].length === 0 ? "": "-L") $(" -L ".join(inputs["intervalList"].map(function(el) { return el.path}))) \
            $(inputs["excludeIntervalList"].length === 0 ? "": "-XL") $(" -XL ".join(inputs["excludeIntervalList"].map(function(el) { return el.path}))) \
            $(inputs["dontUseSoftClippedBases"] ? "--dont-use-soft-clipped-bases" : "") \

  - class: InlineJavascriptRequirement
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes* 60)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
