cwlVersion: v1.2
id: blast
class: Workflow
requirements:
  - class: InlineJavascriptRequirement
inputs:
  - id: blast_docker_override
    type:
      - string
      - 'null'
  - id: blast_docker
    type:
      - string
      - 'null'
  - id: queryfa
    type: File
  - id: fname
    default: /sfs/blastdb/2019_ncov/nucl/v6/ncov
    type: string
  - id: method
    default: blastn
    type: string
  - id: outfmt
    default: 7
    type: int
  - id: evalue
    default: 10.0
    type: float
  - id: Outfile
    type:
      - string
      - 'null'
  - id: threads
    default: 8
    type: int
  - id: runblastp.max_target_seqs
    default: 100
    type: int
  - id: runblastp.word_size
    default: 6
    type: int
  - id: runblastp.seg
    default: no
    type: string
  - id: runblastp.comp_based_stats
    default: '2'
    type: string
  - id: runblastp.matrix
    default: BLOSUM62
    type: string
  - id: runblastp.gapopen
    default: 11
    type: int
  - id: runblastp.gapextend
    default: 1
    type: int
  - id: runblastp.max_hsps
    type:
      - int
      - 'null'
  - id: runblastp.taxids
    type:
      - string
      - 'null'
  - id: runblastp.negative_taxids
    type:
      - string
      - 'null'
  - id: runblastp.lcase_masking
    default: false
    type: boolean
  - id: runblastn.max_target_seqs
    default: 100
    type: int
  - id: runblastn.word_size
    default: 28
    type: int
  - id: runblastn.reward
    default: 1
    type: int
  - id: runblastn.penalty
    default: -2
    type: int
  - id: runblastn.strand
    default: both
    type: string
  - id: runblastn.gapopen
    default: 0
    type: int
  - id: runblastn.gapextend
    default: 0
    type: int
  - id: runblastn.dust
    default: "'20 64 1'"
    type: string
  - id: runblastn.max_hsps
    type:
      - int
      - 'null'
  - id: runblastn.tasks
    default: megablast
    type: string
  - id: runblastn.taxids
    type:
      - string
      - 'null'
  - id: runblastn.negative_taxids
    type:
      - string
      - 'null'
  - id: runblastn.lcase_masking
    default: false
    type: boolean
  - id: runblastx.max_target_seqs
    default: 100
    type: int
  - id: runblastx.word_size
    default: 6
    type: int
  - id: runblastx.seg
    default: "'12 2.2 2.5'"
    type: string
  - id: runblastx.comp_based_stats
    default: '2'
    type: string
  - id: runblastx.matrix
    default: BLOSUM62
    type: string
  - id: runblastx.gapopen
    default: 11
    type: int
  - id: runblastx.gapextend
    default: 1
    type: int
  - id: runblastx.taxids
    type:
      - string
      - 'null'
  - id: runblastx.negative_taxids
    type:
      - string
      - 'null'
  - id: runblastx.max_hsps
    type:
      - int
      - 'null'
  - id: runblastx.lcase_masking
    default: false
    type: boolean
  - id: runtblastn.max_target_seqs
    default: 100
    type: int
  - id: runtblastn.word_size
    default: 6
    type: int
  - id: runtblastn.seg
    default: "'12 2.2 2.5'"
    type: string
  - id: runtblastn.comp_based_stats
    default: '2'
    type: string
  - id: runtblastn.matrix
    default: BLOSUM62
    type: string
  - id: runtblastn.gapopen
    default: 11
    type: int
  - id: runtblastn.gapextend
    default: 1
    type: int
  - id: runtblastn.lcase_masking
    default: false
    type: boolean
  - id: runtblastn.max_hsps
    type:
      - int
      - 'null'
  - id: runtblastn.taxids
    type:
      - string
      - 'null'
  - id: runtblastn.negative_taxids
    type:
      - string
      - 'null'
  - id: runtblastx.taxids
    type:
      - string
      - 'null'
  - id: runtblastx.word_size
    default: 3
    type: int
  - id: runtblastx.max_target_seqs
    default: 100
    type: int
  - id: runtblastx.seg
    default: "'12 2.2 2.5'"
    type: string
  - id: runtblastx.matrix
    default: BLOSUM62
    type: string
  - id: runtblastx.lcase_masking
    default: false
    type: boolean
  - id: runtblastx.negative_taxids
    type:
      - string
      - 'null'
  - id: runtblastx.max_hsps
    type:
      - int
      - 'null'
steps:
  - id: _fina_output_select_first
    in:
      - id: runtblastx_out
        source: runtblastx/out
      - id: runblastp_out
        source: runblastp/out
      - id: runblastn_out
        source: runblastn/out
      - id: runblastx_out
        source: runblastx/out
      - id: runtblastn_out
        source: runtblastn/out
    out:
      - result
    run:
        class: ExpressionTool
        inputs:
          - id: runtblastx_out
            type: Any
          - id: runblastp_out
            type: Any
          - id: runblastn_out
            type: Any
          - id: runblastx_out
            type: Any
          - id: runtblastn_out
            type: Any
        outputs:
          - id: result
            type: File
        expression: '${ return {"result": [runtblastx_out, runblastp_out, runblastn_out,
            runblastx_out, runtblastn_out].find(function(element) { return element
            !== null }) }; }'
  - id: runblastp
    in:
      - id: docker
        source: blast_docker
      - id: Queryfa
        source: queryfa
      - id: Fname
        source: fname
      - id: outfmt
        source: outfmt
      - id: Outfile
        source: Outfile
      - id: evalue
        source: evalue
      - id: threads
        source: threads
      - id: max_target_seqs
        source: runblastp.max_target_seqs
      - id: word_size
        source: runblastp.word_size
      - id: seg
        source: runblastp.seg
      - id: comp_based_stats
        source: runblastp.comp_based_stats
      - id: matrix
        source: runblastp.matrix
      - id: gapopen
        source: runblastp.gapopen
      - id: gapextend
        source: runblastp.gapextend
      - id: max_hsps
        source: runblastp.max_hsps
      - id: taxids
        source: runblastp.taxids
      - id: negative_taxids
        source: runblastp.negative_taxids
      - id: lcase_masking
        source: runblastp.lcase_masking
    out:
      - id: out
    run:
        class: CommandLineTool
        id: runblastp
        inputs:
          - id: docker
            type: string
          - id: Queryfa
            type: File
          - id: Fname
            type: string
          - id: outfmt
            type: int
          - id: Outfile
            type: string
          - id: evalue
            type: float
          - id: threads
            type: int
          - id: max_target_seqs
            default: 100
            type: int
          - id: word_size
            default: 6
            type: int
          - id: seg
            default: no
            type: string
          - id: comp_based_stats
            default: '2'
            type: string
          - id: matrix
            default: BLOSUM62
            type: string
          - id: gapopen
            default: 11
            type: int
          - id: gapextend
            default: 1
            type: int
          - id: max_hsps
            type:
              - int
              - 'null'
          - id: taxids
            type:
              - string
              - 'null'
          - id: negative_taxids
            type:
              - string
              - 'null'
          - id: lcase_masking
            default: false
            type: boolean
        outputs:
          - id: out
            type: File
            outputBinding:
                glob: $(inputs.Outfile)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4+

                    set -e
                    blastp -db "$(inputs.Fname)" \
                    -query $(inputs.Queryfa.path) \
                    -outfmt $(inputs.outfmt) \
                    -out     $(inputs.Outfile) \
                    -max_target_seqs $(inputs.max_target_seqs) \
                    -comp_based_stats $(inputs.comp_based_stats) \
                    -evalue $(inputs.evalue) \
                    -word_size $(inputs.word_size) \
                    -matrix   $(inputs.matrix) \
                    -seg     $(inputs.seg) \
                    -gapopen $(inputs.gapopen) \
                    -gapextend $(inputs.gapextend) \
                    -num_threads $(inputs.threads) \
                    $(inputs.lcase_masking ? "-lcase_masking" : "") $(inputs.max_hsps === null ? "" : "-max_hsps " + inputs.max_hsps) $(inputs.taxids === null ? "" : "-taxids " + inputs.taxids) $(inputs.negative_taxids === null ? "" : "-negative_taxids " + inputs.negative_taxids) \

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 15258.7890625
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
    when: $(inputs.method === "blastp")
  - id: runblastn
    in:
      - id: docker
        source: blast_docker
      - id: Queryfa
        source: queryfa
      - id: Fname
        source: fname
      - id: Outfile
        source: Outfile
      - id: threads
        source: threads
      - id: outfmt
        source: outfmt
      - id: max_target_seqs
        source: runblastn.max_target_seqs
      - id: evalue
        source: evalue
      - id: word_size
        source: runblastn.word_size
      - id: reward
        source: runblastn.reward
      - id: penalty
        source: runblastn.penalty
      - id: strand
        source: runblastn.strand
      - id: gapopen
        source: runblastn.gapopen
      - id: gapextend
        source: runblastn.gapextend
      - id: dust
        source: runblastn.dust
      - id: max_hsps
        source: runblastn.max_hsps
      - id: tasks
        source: runblastn.tasks
      - id: taxids
        source: runblastn.taxids
      - id: negative_taxids
        source: runblastn.negative_taxids
      - id: lcase_masking
        source: runblastn.lcase_masking
    out:
      - id: out
    run:
        class: CommandLineTool
        id: runblastn
        inputs:
          - id: docker
            type: string
          - id: Queryfa
            type: File
          - id: Fname
            type: string
          - id: Outfile
            type: string
          - id: threads
            type: int
          - id: outfmt
            type: int
          - id: max_target_seqs
            default: 100
            type: int
          - id: evalue
            type: float
          - id: word_size
            default: 28
            type: int
          - id: reward
            default: 1
            type: int
          - id: penalty
            default: -2
            type: int
          - id: strand
            default: both
            type: string
          - id: gapopen
            default: 0
            type: int
          - id: gapextend
            default: 0
            type: int
          - id: dust
            default: "'20 64 1'"
            type: string
          - id: max_hsps
            type:
              - int
              - 'null'
          - id: tasks
            default: megablast
            type: string
          - id: taxids
            type:
              - string
              - 'null'
          - id: negative_taxids
            type:
              - string
              - 'null'
          - id: lcase_masking
            default: false
            type: boolean
        outputs:
          - id: out
            type: File
            outputBinding:
                glob: $(inputs.Outfile)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4+

                    set -e
                    blastn -db "$(inputs.Fname)" \
                    -show_gis \
                    -query $(inputs.Queryfa.path) \
                    -outfmt $(inputs.outfmt) \
                    -out     $(inputs.Outfile) \
                    -max_target_seqs $(inputs.max_target_seqs) \
                    -evalue $(inputs.evalue) \
                    -word_size $(inputs.word_size) \
                    -penalty $(inputs.penalty) \
                    -reward  $(inputs.reward) \
                    -dust $(inputs.dust) \
                    -gapopen $(inputs.gapopen) \
                    -gapextend $(inputs.gapextend) \
                    -task $(inputs.tasks) \
                    -strand  $(inputs.strand) \
                    -num_threads $(inputs.threads) \
                    $(inputs.lcase_masking ? "-lcase_masking" : "") $(inputs.max_hsps === null ? "" : "-max_hsps " + inputs.max_hsps) $(inputs.taxids === null ? "" : "-taxids " + inputs.taxids) $(inputs.negative_taxids === null ? "" : "-negative_taxids " + inputs.negative_taxids)\

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 15258.7890625
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
    when: $(inputs.method === "blastn")
  - id: runblastx
    in:
      - id: Queryfa
        source: queryfa
      - id: Fname
        source: fname
      - id: outfmt
        source: outfmt
      - id: evalue
        source: evalue
      - id: Outfile
        source: Outfile
      - id: docker
        source: blast_docker
      - id: threads
        source: threads
      - id: max_target_seqs
        source: runblastx.max_target_seqs
      - id: word_size
        source: runblastx.word_size
      - id: seg
        source: runblastx.seg
      - id: comp_based_stats
        source: runblastx.comp_based_stats
      - id: matrix
        source: runblastx.matrix
      - id: gapopen
        source: runblastx.gapopen
      - id: gapextend
        source: runblastx.gapextend
      - id: taxids
        source: runblastx.taxids
      - id: negative_taxids
        source: runblastx.negative_taxids
      - id: max_hsps
        source: runblastx.max_hsps
      - id: lcase_masking
        source: runblastx.lcase_masking
    out:
      - id: out
    run:
        class: CommandLineTool
        id: runblastx
        inputs:
          - id: Queryfa
            type: File
          - id: Fname
            type: string
          - id: outfmt
            type: int
          - id: evalue
            type: float
          - id: Outfile
            type: string
          - id: docker
            type: string
          - id: threads
            type: int
          - id: max_target_seqs
            default: 100
            type: int
          - id: word_size
            default: 6
            type: int
          - id: seg
            default: "'12 2.2 2.5'"
            type: string
          - id: comp_based_stats
            default: '2'
            type: string
          - id: matrix
            default: BLOSUM62
            type: string
          - id: gapopen
            default: 11
            type: int
          - id: gapextend
            default: 1
            type: int
          - id: taxids
            type:
              - string
              - 'null'
          - id: negative_taxids
            type:
              - string
              - 'null'
          - id: max_hsps
            type:
              - int
              - 'null'
          - id: lcase_masking
            default: false
            type: boolean
        outputs:
          - id: out
            type: File
            outputBinding:
                glob: $("$" + inputs.Outfile)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4+

                    set -e
                    blastx -db "$(inputs.Fname)" \
                    -query $(inputs.Queryfa.path) \
                    -outfmt $(inputs.outfmt) \
                    -out     $(inputs.Outfile) \
                    -max_target_seqs $(inputs.max_target_seqs) \
                    -comp_based_stats $(inputs.comp_based_stats) \
                    -evalue $(inputs.evalue) \
                    -word_size $(inputs.word_size) \
                    -matrix   $(inputs.matrix) \
                    -seg     $(inputs.seg) \
                    -gapopen $(inputs.gapopen) \
                    -gapextend $(inputs.gapextend) \
                    -num_threads $(inputs.threads) \
                    $(inputs.lcase_masking ? "-lcase_masking" : "") $(inputs.max_hsps === null ? "" : "-max_hsps " + inputs.max_hsps) $(inputs.taxids === null ? "" : "-taxids " + inputs.taxids) $(inputs.negative_taxids === null ? "" : "-negative_taxids " + inputs.negative_taxids)\

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 15258.7890625
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
    when: $(inputs.method === "blastx")
  - id: runtblastn
    in:
      - id: Queryfa
        source: queryfa
      - id: Fname
        source: fname
      - id: outfmt
        source: outfmt
      - id: evalue
        source: evalue
      - id: Outfile
        source: Outfile
      - id: docker
        source: blast_docker
      - id: threads
        source: threads
      - id: max_target_seqs
        source: runtblastn.max_target_seqs
      - id: word_size
        source: runtblastn.word_size
      - id: seg
        source: runtblastn.seg
      - id: comp_based_stats
        source: runtblastn.comp_based_stats
      - id: matrix
        source: runtblastn.matrix
      - id: gapopen
        source: runtblastn.gapopen
      - id: gapextend
        source: runtblastn.gapextend
      - id: lcase_masking
        source: runtblastn.lcase_masking
      - id: max_hsps
        source: runtblastn.max_hsps
      - id: taxids
        source: runtblastn.taxids
      - id: negative_taxids
        source: runtblastn.negative_taxids
    out:
      - id: out
    run:
        class: CommandLineTool
        id: runtblastn
        inputs:
          - id: Queryfa
            type: File
          - id: Fname
            type: string
          - id: outfmt
            type: int
          - id: evalue
            type: float
          - id: Outfile
            type: string
          - id: docker
            type: string
          - id: threads
            type: int
          - id: max_target_seqs
            default: 100
            type: int
          - id: word_size
            default: 6
            type: int
          - id: seg
            default: "'12 2.2 2.5'"
            type: string
          - id: comp_based_stats
            default: '2'
            type: string
          - id: matrix
            default: BLOSUM62
            type: string
          - id: gapopen
            default: 11
            type: int
          - id: gapextend
            default: 1
            type: int
          - id: lcase_masking
            default: false
            type: boolean
          - id: max_hsps
            type:
              - int
              - 'null'
          - id: taxids
            type:
              - string
              - 'null'
          - id: negative_taxids
            type:
              - string
              - 'null'
        outputs:
          - id: out
            type: File
            outputBinding:
                glob: $(inputs.Outfile)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4+

                    set -e
                    tblastn -db "$(inputs.Fname)" \
                    -query $(inputs.Queryfa.path) \
                    -outfmt $(inputs.outfmt) \
                    -out    $(inputs.Outfile) \
                    -max_target_seqs $(inputs.max_target_seqs) \
                    -comp_based_stats $(inputs.comp_based_stats) \
                    -evalue $(inputs.evalue) \
                    -word_size $(inputs.word_size) \
                    -matrix   $(inputs.matrix) \
                    -seg     $(inputs.seg) \
                    -gapopen $(inputs.gapopen) \
                    -gapextend $(inputs.gapextend) \
                    -num_threads $(inputs.threads) \
                    $(inputs.lcase_masking ? "-lcase_masking" : "") $(inputs.max_hsps === null ? "" : "-max_hsps " + inputs.max_hsps) $(inputs.taxids === null ? "" : "-taxids " + inputs.taxids) $(inputs.negative_taxids === null ? "" : "-negative_taxids " + inputs.negative_taxids)\

          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 15258.7890625
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
    when: $(inputs.method === "queryfa")
  - id: runtblastx
    in:
      - id: Queryfa
        source: queryfa
      - id: Fname
        source: fname
      - id: outfmt
        source: outfmt
      - id: Outfile
        source: Outfile
      - id: threads
        source: threads
      - id: evalue
        source: evalue
      - id: docker
        source: blast_docker
      - id: taxids
        source: runtblastx.taxids
      - id: word_size
        source: runtblastx.word_size
      - id: max_target_seqs
        source: runtblastx.max_target_seqs
      - id: seg
        source: runtblastx.seg
      - id: matrix
        source: runtblastx.matrix
      - id: lcase_masking
        source: runtblastx.lcase_masking
      - id: negative_taxids
        source: runtblastx.negative_taxids
      - id: max_hsps
        source: runtblastx.max_hsps
    out:
      - id: out
    run:
        class: CommandLineTool
        id: runtblastx
        inputs:
          - id: Queryfa
            type: File
          - id: Fname
            type: string
          - id: outfmt
            type: int
          - id: Outfile
            type: string
          - id: threads
            type: int
          - id: evalue
            type: float
          - id: docker
            type: string
          - id: taxids
            type:
              - string
              - 'null'
          - id: word_size
            default: 3
            type: int
          - id: max_target_seqs
            default: 100
            type: int
          - id: seg
            default: "'12 2.2 2.5'"
            type: string
          - id: matrix
            default: BLOSUM62
            type: string
          - id: lcase_masking
            default: false
            type: boolean
          - id: negative_taxids
            type:
              - string
              - 'null'
          - id: max_hsps
            type:
              - int
              - 'null'
        outputs:
          - id: out
            type: File
            outputBinding:
                glob: $(inputs.Outfile)
        requirements:
          - class: InitialWorkDirRequirement
            listing:
              - entryname: script.bash
                entry: |4

                    set -e
                    tblastx -db "$(inputs.Fname)" \
                    -query $(inputs.Queryfa.path) \
                    -outfmt $(inputs.outfmt) \
                    -out     $(inputs.Outfile) \
                    -max_target_seqs $(inputs.max_target_seqs) \
                    -evalue $(inputs.evalue) \
                    -word_size $(inputs.word_size) \
                    -matrix   $(inputs.matrix) \
                    -seg     $(inputs.seg) \
                    -num_threads $(inputs.threads) \
                    $(inputs.lcase_masking ? "-lcase_masking" : "") $(inputs.max_hsps === null ? "" : "-max_hsps " + inputs.max_hsps) $(inputs.taxids === null ? "" : "-taxids " + inputs.taxids) $(inputs.negative_taxids === null ? "" : "-negative_taxids " + inputs.negative_taxids)\
          - class: InlineJavascriptRequirement
          - class: NetworkAccess
            networkAccess: true
        hints:
          - class: ResourceRequirement
            coresMin: 8
            ramMin: 15258.7890625
            outdirMin: 1024
        cwlVersion: v1.2
        baseCommand:
          - bash
          - script.bash
    when: $(inputs.method === "tblastx")
outputs:
  - id: blast.fina_output
    outputSource: _fina_output_select_first/result
    type: File
