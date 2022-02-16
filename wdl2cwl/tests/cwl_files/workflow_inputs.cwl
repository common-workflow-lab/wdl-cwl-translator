cwlVersion: v1.2
id: foo
class: Workflow
inputs:
  - id: first
    doc: test coverage example
    default: 'one '
    type: string
  - id: third
    default: 3.14159
    type: float
steps: []
outputs:
  - id: foo.first_result
    outputSource: first
    type: string
  - id: foo.third_result
    outputSource: third
    type: float
