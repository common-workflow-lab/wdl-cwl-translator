cwlVersion: v1.2
id: quoteWorkflow
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
inputs:
  - id: str_arr
    type:
      name: _str_arr_string_array
      items: string
      type: array
  - id: int_arr
    type:
      name: _int_arr_int_array
      items: int
      type: array
  - id: float_arr
    type:
      name: _float_arr_float_array
      items: float
      type: array
  - id: bool_arr
    type:
      name: _bool_arr_boolean_array
      items: boolean
      type: array
  - id: file_arr
    type:
      name: _file_arr_File_array
      items: File
      type: array
steps:
  - id: _str_output_quote
    in:
      - id: str_arr
        source: str_arr
    out:
      - result
    run:
      id: _str_output_quote_etool
      class: ExpressionTool
      inputs:
        - id: str_arr
          type: Any
      outputs:
        - id: result
          type:
            name: _str_output_string_array
            items: string
            type: array
      expression: "${ return {\"result\": inputs.str_arr.map(function(item) {return
        '\\\"'+item+'\\\"'})}; }"
  - id: _int_output_quote
    in:
      - id: int_arr
        source: int_arr
    out:
      - result
    run:
      id: _int_output_quote_etool
      class: ExpressionTool
      inputs:
        - id: int_arr
          type: Any
      outputs:
        - id: result
          type:
            name: _int_output_string_array
            items: string
            type: array
      expression: "${ return {\"result\": inputs.int_arr.map(function(item) {return
        '\\\"'+item+'\\\"'})}; }"
  - id: _float_output_quote
    in:
      - id: float_arr
        source: float_arr
    out:
      - result
    run:
      id: _float_output_quote_etool
      class: ExpressionTool
      inputs:
        - id: float_arr
          type: Any
      outputs:
        - id: result
          type:
            name: _float_output_string_array
            items: string
            type: array
      expression: "${ return {\"result\": inputs.float_arr.map(function(item) {return
        '\\\"'+item+'\\\"'})}; }"
  - id: _bool_output_quote
    in:
      - id: bool_arr
        source: bool_arr
    out:
      - result
    run:
      id: _bool_output_quote_etool
      class: ExpressionTool
      inputs:
        - id: bool_arr
          type: Any
      outputs:
        - id: result
          type:
            name: _bool_output_string_array
            items: string
            type: array
      expression: "${ return {\"result\": inputs.bool_arr.map(function(item) {return
        '\\\"'+item+'\\\"'})}; }"
  - id: file_quote
    in:
      - id: file_arr
        source: file_arr
    out:
      - id: out
    run:
      id: file_quote
      class: CommandLineTool
      inputs:
        - id: file_arr
          type:
            name: _file_arr_File_array
            items: File
            type: array
      outputs:
        - id: out
          type: File
          outputBinding:
            glob: output.txt
      requirements:
        - class: InitialWorkDirRequirement
          listing:
            - entryname: script.bash
              entry: |2

                echo $(inputs.file_arr.map(function(item) {return '\"'+item.path+'\"'}).join(" ")) >> output.txt
        - class: InlineJavascriptRequirement
        - class: NetworkAccess
          networkAccess: true
      hints:
        - class: ResourceRequirement
          outdirMin: 1024
      cwlVersion: v1.2
      baseCommand:
        - bash
        - script.bash
outputs:
  - id: quoteWorkflow.str_output
    outputSource: _str_output_quote/result
    type:
      name: _str_output_string_array
      items: string
      type: array
  - id: quoteWorkflow.int_output
    outputSource: _int_output_quote/result
    type:
      name: _int_output_string_array
      items: string
      type: array
  - id: quoteWorkflow.float_output
    outputSource: _float_output_quote/result
    type:
      name: _float_output_string_array
      items: string
      type: array
  - id: quoteWorkflow.bool_output
    outputSource: _bool_output_quote/result
    type:
      name: _bool_output_string_array
      items: string
      type: array
  - id: quoteWorkflow.file_output
    outputSource: file_quote/out
    type: File
