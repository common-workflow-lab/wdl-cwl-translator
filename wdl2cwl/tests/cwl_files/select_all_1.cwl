class: CommandLineTool
id: test_select_all
inputs:
  - id: one
    type: int
  - id: two
    type:
      - int
      - 'null'
outputs:
  - id: first1
    type:
        items: int
        type: array
    outputBinding:
        outputEval: $([inputs.one, inputs.two, 1, [inputs.two, inputs.one].find(function(element)
            { return element !== null }) ].filter(function(element) { return element
            !== null }) )
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: ''
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    outdirMin: 1024
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
