cwlVersion: v1.2
id: foo
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
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
        valueFrom: $(self.one + ".suffix")
      - id: other
        source: fourth
        valueFrom: $(self.two)
      - id: echo
        default: true
    out:
      - id: out
      - id: result
    run:
        id: echo
        class: CommandLineTool
        inputs:
          - id: in
            type: string
          - id: other
            type: string
          - id: echo
            type: boolean
        outputs:
          - id: out
            type: string
            outputBinding:
                outputEval: $(inputs["in"])
          - id: result
            type: string
            outputBinding:
                loadContents: true
                glob: _stdout
                outputEval: $(self[0].contents.replace(/[\r\n]+$/, ''))
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    $(inputs.echo ? "echo " + inputs["in"] + " " + inputs.other : "")
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
        stdout: _stdout
outputs:
  - id: foo.first_result
    outputSource: first
    type: string
  - id: foo.third_result
    outputSource: third
    type: float
  - id: foo.echo_out
    outputSource: echo/out
    type: string
  - id: foo.echo_result
    outputSource: echo/result
    type: string
