cwlVersion: v1.2
$graph:
- class: Workflow
  id: main
  inputs:
    - id: blast_docker_override
      type:
        - string
        - 'null'
    - id: blast_docker
      type: string
      # TODO: how to translate the original?
      #       String blast_docker = select_first([blast_docker_override,"swr.cn-south-1.myhuaweicloud.com/cngbdb/blast:1.2"])
    - id: queryfa
      type: File
    - id: fname
      type: string
      default: '/sfs/blastdb/2019_ncov/nucl/v6/ncov'
    - id: method
      type: string
      default: 'blastn'
    - id: outfmt
      type: int
      default: 7
    - id: evalue
      type: float
      default: 10
    # TODO: how to translate the original?
    #       String Outfile = basename(queryfa)+'.blast_result.txt'
    - id: Outfile
      type: string
    # runblastn task inputs, to be linked from workflow -> tool
    - id: runblastn.docker
      type: string
    - id: runblastn.Queryfa
      type: File
    - id: runblastn.Fname
      type: string
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
    # runblastp task inputs, to be linked from workflow -> tool
    - id: runblastp.docker
      type: string
    - id: runblastp.Queryfa
      type: File
    - id: runblastp.Fname
      type: string
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
  steps:
    - id: runblastp
      in:
        Fname: fname
        Queryfa: queryfa
        docker: blast_docker
        outfmt: outfmt
        Outfile: Outfile
        threads: threads
      run: runblastp
      when: $(inputs.method == 'blastp')
      out: [out]
    - id: runblastn
      in:
        Fname: fname
        Queryfa: queryfa
        docker: blast_docker
        outfmt: outfmt
        evalue: evalue
        Outfile: Outfile
        threads: threads
      run: runblastn
      when: $(inputs.method == 'blastn')
      out: [out]
  outputs:
    - id: fina_output
      type: File
      outputSource:
        - runblastp/out
        - runblastn/out
      pickValue: first_non_null

- class: CommandLineTool
  id: runblastn
  requirements:
    - class: InitialWorkDirRequirement
      listing:
        - entryname: script.bash
          entry: |4
              set -e
              blastn -db "$(inputs.Fname)" \
              -show_gis \
              -query $(inputs.QueryFa) \
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
              $(true='-lcase_masking' false='' inputs.lcase_masking) $("-max_hsps "+ inputs.max_hsps) $("-taxids " +inputs.taxids) $("-negative_taxids " +inputs.negative_taxids)\
  hints:
    - class: DockerRequirement
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 16384
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
  baseCommand:
    - bash
    - script.bash
  outputs:
    - id: out
      type: File
      outputBinding:
        glob: $(inputs.Outfile)

- class: CommandLineTool
  id: runblastp
  requirements:
    - class: InitialWorkDirRequirement
      listing:
        - entryname: script.bash
          entry: |4
              set -e
              blastp -db "$(inputs.Fname)" \
              -query $(inputs.Queryfa) \
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
              $(true='-lcase_masking' false='' inputs.lcase_masking) $("-max_hsps "+inputs.max_hsps) $("-taxids " +inputs.taxids) $("-negative_taxids " +inputs.negative_taxids) \
  hints:
    - class: DockerRequirement
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 16384
  inputs:
    - id: docker
      type: string
    - id: Queryfa
      type: File
    - id: Fname
      type: string
    - id: outfmt
      type: string
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
  baseCommand:
    - bash
    - script.bash
  outputs:
    - id: out
      type: File
      outputBinding:
        glob: $(inputs.Outfile)
