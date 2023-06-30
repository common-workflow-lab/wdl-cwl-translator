cwlVersion: v1.2
id: literal
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            echo $(inputs.tenth.map(function(el) {return el.path}).join(" "))
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: ResourceRequirement
    outdirMin: 1024
inputs:
  - id: first
    default: true
    type: boolean
  - id: second
    default: 42
    type: int
  - id: third
    default: 6.022
    type: float
  - id: fourth
    default: hoopla
    type: string
  - id: fifth
    default:
        class: File
        path: ../../../README.md
    type: File
  - id: sixth
    default:
      - true
      - false
    type:
        items: boolean
        type: array
  - id: seventh
    default:
      - 42
      - 23
    type:
        items: int
        type: array
  - id: eighth
    default:
      - 6.022
      - 10.0
      - 23.0
    type:
        items: float
        type: array
  - id: nineth
    default:
      - Hello
      - World
    type:
        items: string
        type: array
  - id: tenth
    default:
      - class: File
        path: ../../../README.md
      - class: File
        path: ../../../LICENSE
    type:
        items: File
        type: array
baseCommand:
  - bash
  - script.bash
outputs:
  - id: result
    type: string
    outputBinding:
        outputEval: $(inputs.first + " " + inputs.second + " " + inputs.third + "
            " + inputs.fourth + " " + inputs.fifth.basename + " " + inputs.sixth.join(",")
            + " " + inputs.seventh.join(",") + " " + inputs.eighth.join(",") + " "
            + inputs.nineth.join(",") + " " + inputs.tenth[0].basename + "," + inputs.tenth[1].basename)
