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
      name: _a_s_string_array
      items: string
      type: array
  - id: a_s2
    type:
      name: _a_s2_string_array
      items: string
      type: array
baseCommand:
  - 'true'
outputs:
  - id: out_s
    type:
      name: _out_s_string_array
      items: string
      type: array
    outputBinding:
      outputEval: $((function () {var new_array = []; [ inputs.a_s, inputs.a_s2 
        ].forEach(function(value, index, obj) {value.forEach(function(sub_value,
        sub_index, sub_obj) {new_array.push(sub_value);});}); return 
        new_array;})())
  - id: out_a_s
    type:
      name: _out_a_s_string_array
      items: string
      type: array
    outputBinding:
      outputEval: $(inputs.a_s)
