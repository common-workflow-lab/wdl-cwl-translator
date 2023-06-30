cwlVersion: v1.2
id: literal_struct_test
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |2

            echo $({ "one": "four", "two": "five" }.one) $({ "one": "four", "two": "five" }.two) $(inputs.six.one) $(inputs.six.two)
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: ResourceRequirement
    outdirMin: 1024
inputs:
  - id: six
    default:
        one: seven
        two: eight
    type:
        name: Foo
        fields:
          - name: one
            type: string
          - name: two
            type: string
        type: record
baseCommand:
  - bash
  - script.bash
stdout: _stdout
outputs:
  - id: result
    type: string
    outputBinding:
        loadContents: true
        glob: _stdout
        outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
