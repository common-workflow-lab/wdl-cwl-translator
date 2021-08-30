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
            gatk --java-options '-Xmx$(inputs.javaXmxMb)M -XX:ParallelGCThreads=1' \
            HaplotypeCaller \
            -R wd2/$(inputs.basename(referenceFasta)) \
            -O $(inputs.outputPath) \
            -I ${
            var text = "";
            for(var i=0;i<inputs["inputBams"].length;i++) 
              text+= inputs["inputBams"][i].path+" -I ";
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
