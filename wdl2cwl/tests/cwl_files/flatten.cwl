class: CommandLineTool
id: echo
inputs:
  - id: a_s
    type:
        items: string
        type: array
  - id: a_s2
    type:
        items: string
        type: array
outputs:
  - id: out_s
    type:
        items: string
        type: array
    outputBinding:
        outputEval: $((function () {var new_array = []; [ inputs.a_s, inputs.a_s2
            ].forEach(function(value, index, obj) {value.forEach(function(sub_value,
            sub_index, sub_obj) {new_array.push(sub_value);});}); return new_array;})())
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.bash
        entry: |4+

  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
  - class: ResourceRequirement
    outdirMin: 1024
cwlVersion: v1.2
baseCommand:
  - bash
  - script.bash
