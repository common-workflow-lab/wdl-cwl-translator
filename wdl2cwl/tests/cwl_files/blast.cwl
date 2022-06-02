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
steps: []
outputs:
  - id: blast.fina_output
    outputSource: '$([inputs["runtblastx.out"] === null '
    type: File
