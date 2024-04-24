cwlVersion: v1.2
id: sepWorkflow
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
inputs:
  - id: delimiter
    type: string
  - id: to_sep
    type:
        items: string
        type: array
steps:
  - id: _str_output_sep
    in:
      - id: delimiter
        source: delimiter
      - id: to_sep
        source: to_sep
    out:
      - result
    run:
        class: ExpressionTool
        inputs:
          - id: delimiter
            type: Any
          - id: to_sep
            type: Any
        outputs:
          - id: result
            type: string
        expression: '${ return {"result": inputs.to_sep.join(inputs.delimiter)}; }'
outputs:
  - id: sepWorkflow.str_output
    outputSource: _str_output_sep/result
    type: string
