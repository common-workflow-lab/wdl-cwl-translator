cwlVersion: v1.2
id: test
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
inputs:
  - id: gtf_version
    type: string
  - id: organism
    type: string
  - id: organism_prefix
    type: string
steps:
  - id: _GetReferences.references.genome_fa
    in:
      - id: target
        source: GetReferences/references
        valueFrom: $(self.genome_fa)
    out:
      - genome_fa
    run:
        class: ExpressionTool
        inputs:
          - id: target
            type: Any
        outputs:
          - id: genome_fa
            type: File
        expression: '${return {"genome_fa": inputs.target}; }'
  - id: _GetReferences.references.annotation_gtf
    in:
      - id: target
        source: GetReferences/references
        valueFrom: $(self.annotation_gtf)
    out:
      - annotation_gtf
    run:
        class: ExpressionTool
        inputs:
          - id: target
            type: Any
        outputs:
          - id: annotation_gtf
            type: File
        expression: '${return {"annotation_gtf": inputs.target}; }'
  - id: GetReferences
    in:
      - id: gtf_version
        source: gtf_version
      - id: organism
        source: organism
      - id: organism_prefix
        source: organism_prefix
    out:
      - id: references
    run:
        class: CommandLineTool
        id: GetReferences
        doc: Download files needed for building the designated references
        inputs:
          - id: gtf_version
            type: string
          - id: organism
            type: string
          - id: organism_prefix
            type: string
        outputs:
          - id: references
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
            outputBinding:
                outputEval: '$({ "genome_fa": { "class": "File", "path": runtime.outdir+"/"+"GRC"
                    + inputs.organism_prefix + "38.primary_assembly.genome.fa" },
                    "annotation_gtf": { "class": "File", "path": runtime.outdir+"/"+"gencode.v"
                    + inputs.gtf_version + ".primary_assembly.annotation.gtf" } })'
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -eo pipefail

                    echo a > $("GRC" + inputs.organism_prefix + "38.primary_assembly.genome.fa")
                    echo b > $("gencode.v" + inputs.gtf_version + ".primary_assembly.annotation.gtf")
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
  - id: EchoRef
    in:
      - id: ref
        source: GetReferences/references
    out:
      - id: refflat
    run:
        class: CommandLineTool
        id: EchoRef
        inputs:
          - id: ref
            type:
                name: References
                fields:
                  - name: genome_fa
                    type: File
                  - name: annotation_gtf
                    type: File
                type: record
        outputs:
          - id: refflat
            type: File
            outputBinding:
                glob: $(inputs.ref.annotation_gtf.basename.replace(/\.gtf$/, '')  +
                    ".refflat.txt")
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    touch $(inputs.ref.annotation_gtf.basename.replace(/\.gtf$/, '')  + ".refflat.txt")
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
  - id: test.genome_fa
    outputSource: _GetReferences.references.genome_fa/genome_fa
    type: File
  - id: test.annotation_gtf
    outputSource: _GetReferences.references.annotation_gtf/annotation_gtf
    type: File
  - id: test.refflat
    outputSource: EchoRef/refflat
    type: File
