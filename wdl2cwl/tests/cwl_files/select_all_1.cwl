cwlVersion: v1.2
id: test_select_all
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: ResourceRequirement
    outdirMin: 1024
inputs:
  - id: one
    type: int
  - id: two
    type:
      - int
      - 'null'
baseCommand:
  - 'true'
outputs:
  - id: first1
    doc: amalgamation
    type:
        items: int
        type: array
    outputBinding:
        outputEval: $([inputs.one, inputs.two, 1, [inputs.two, inputs.one].find(function(element)
            { return element !== null }) ].filter(function(element) { return element
            !== null }) )
