cwlVersion: v1.2
id: placeholder_options_test
class: CommandLineTool
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4

            echo $(true ? "true single quote: 'foo' " : "error")
            echo $(true ? 'true double quote: "foo" ' : "error")
            echo $(true ? "true mixed quotes: \"foo\" 'bar' " : "error")
            echo $(true ? "true mixed quotes: \"foo\" 'bar' " : "error")
            echo $(inputs.absent === null ? "default single 'quote'" : inputs.absent)
            echo $(inputs.absent === null ? 'default double "quote"' : inputs.absent)
            echo $(inputs.absent === null ? "default \"mixed\" 'quotes'" : inputs.absent)
            echo $(inputs.absent === null ? "default \"mixed\" 'quotes'" : inputs.absent)
            echo $([ 1, 2, 3 ].join(","))
            echo $(inputs.missing_numbers === null ? "success" : inputs.missing_numbers.join(","))
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: ResourceRequirement
    outdirMin: 1024
inputs:
  - id: absent
    type:
      - string
      - 'null'
  - id: missing_numbers
    type:
      - items: int
        type: array
      - 'null'
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
  - id: present
    type: boolean
    outputBinding:
        outputEval: $(inputs.absent !== null)
