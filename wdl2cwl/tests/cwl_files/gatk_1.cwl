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
            for FILE in $(inputs.sep(" ",inputBams)); do ln -s $FILE wd/\$(inputBams $FILE) ; done
            for FILE in $(inputs.sep(" ",inputBamsIndex)); do ln -s $FILE wd/\$(inputBamsIndex $FILE) ; done
            mkdir wd2
            ln -s $(inputs.referenceFasta.path) wd2/$(inputs.basename(referenceFasta))
            ln -s $(inputs.referenceFastaDict.path) wd2/$(inputs.basename(referenceFastaDict))
            ln -s $(inputs.referenceFastaFai) wd2/$(inputs.basename(referenceFastaFai))
            gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
            HaplotypeCaller \
            -R wd2/$(inputs.basename(referenceFasta)) \
            -O $(inputs.outputPath) \
            (for FILE in $(inputs.sep(" ",inputBams)); do echo -- "-I wd/"\$(basename $FILE); done)
            $(inputs["intervalList"].length === 0 ? "": "-L") ${
            var text = "";
            var arr_length = inputs["intervalList"].length;
            for(var i=0;i<arr_length-1;i++) 
              text+= inputs["intervalList"][i].path+" -L ";
            text+= inputs["intervalList"][arr_length-1].path;
            return text;
            } \
            $(inputs["excludeIntervalList"].length === 0 ? "": "-XL") ${
            var text = "";
            var arr_length = inputs["excludeIntervalList"].length;
            for(var i=0;i<arr_length-1;i++) 
              text+= inputs["excludeIntervalList"][i].path+" -XL ";
            text+= inputs["excludeIntervalList"][arr_length-1].path;
            return text;
            } \
            $(inputs["dontUseSoftClippedBases"] ? "--dont-use-soft-clipped-bases" : "") \

  - class: InlineJavascriptRequirement
  - class: ToolTimeLimit
    timelimit: $(inputs.timeMinutes* 60)
cwlVersion: v1.2
baseCommand:
  - sh
  - example.sh
