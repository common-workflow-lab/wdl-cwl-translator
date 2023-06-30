cwlVersion: v1.2
id: nonempty_array_test
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            echo $(inputs.numbers.join(","))
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: ResourceRequirement
    outdirMin: 1024
inputs:
  - id: numbers
    default:
      - 1
      - 2
      - 3
    type:
        items: int
        type: array
baseCommand:
  - bash
  - script.bash
arguments:
  - valueFrom: ${if (inputs.numbers.length == 0) {throw "numbers must contain at least
        one item.";} else { return "";}}
stdout: _stdout
outputs:
  - id: result
    type: string
    outputBinding:
        loadContents: true
        glob: _stdout
        outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
