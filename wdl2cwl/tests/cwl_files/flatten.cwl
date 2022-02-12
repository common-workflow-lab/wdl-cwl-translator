cwlVersion: v1.2
id: echo
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
  - class: NetworkAccess
    networkAccess: true
hints:
  - class: ResourceRequirement
    outdirMin: 1024
inputs:
  - id: a_s
    type:
        items: string
        type: array
  - id: a_s2
    type:
        items: string
        type: array
baseCommand:
  - 'true'
outputs:
  - id: out_s
    type:
        items: string
        type: array
    outputBinding:
        outputEval: $((function () {var new_array = []; [ inputs.a_s, inputs.a_s2
            ].forEach(function(value, index, obj) {value.forEach(function(sub_value,
            sub_index, sub_obj) {new_array.push(sub_value);});}); return new_array;})())
  - id: out_a_s
    type:
        items: string
        type: array
    outputBinding:
        outputEval: $(inputs.a_s)
