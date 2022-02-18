cwlVersion: v1.2
id: foo
class: Workflow
requirements:
  - class: StepInputExpressionRequirement
inputs:
  - id: first
    doc: test coverage example
    default: 'one '
    type: string
  - id: third
    default: 3.14159
    type: float
  - id: fourth
    default:
        one: fifth
        two: sixth
    type:
        name: Foo
        fields:
          - name: one
            type: string
          - name: two
            type: string
        type: record
steps:
  - id: echo
    in:
      - id: in
        source: fourth
        valueFrom: $(self.one)
    out:
      - id: out
    run:
        class: CommandLineTool
        id: echo
        inputs:
          - id: in
            type: string
        outputs:
          - id: out
            type: string
            outputBinding:
                outputEval: $(inputs["in"])
        requirements:
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: ResourceRequirement
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - 'true'
outputs:
  - id: foo.first_result
    outputSource: first
    type: string
  - id: foo.third_result
    outputSource: third
    type: float
  - id: foo.echo_result
    outputSource: echo/out
    type: string
