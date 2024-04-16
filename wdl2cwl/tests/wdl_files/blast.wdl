version 1.0

workflow blast {
  input {
    String? blast_docker_override
    String  blast_docker= select_first([blast_docker_override,"swr.cn-south-1.myhuaweicloud.com/cngbdb/blast:1.2"])
    File    queryfa
    String  fname       = '/sfs/blastdb/2019_ncov/nucl/v6/ncov'
    String  method      = 'blastn'
    Int     outfmt      = 7
    Float   evalue      = 10
    String  Outfile     = basename(queryfa)+'.blast_result.txt'
    Int     threads     = 8
  }
  if (method == 'blastp') {
    call runblastp{
      input:
        Fname      = fname,
        Queryfa    = queryfa,
        docker     = blast_docker,
        outfmt     = outfmt,
        evalue     = evalue,
        Outfile    = Outfile,
        threads    = threads
    }
  }
  if ( method == 'blastn'){
    call runblastn{
      input:
        Fname      = fname,
        Queryfa    = queryfa,
        docker     = blast_docker,
        outfmt     = outfmt,
        evalue     = evalue,
        Outfile    = Outfile,
        threads    = threads
    }
  }
  if ( method == 'blastx'){
    call runblastx{
      input:
        Fname      = fname,
        Queryfa    = queryfa,
        docker     = blast_docker,
        outfmt     = outfmt,
        evalue     = evalue,
        Outfile    = Outfile,
        threads    = threads
    }
  }
  if ( method == 'queryfa'){
    call runtblastn{
      input:
        Fname      = fname,
        Queryfa    = queryfa,
        docker     = blast_docker,
        outfmt     = outfmt,
        evalue     = evalue,
        Outfile    = Outfile,
        threads    = threads
    }
  }
  if ( method == 'tblastx'){
    call runtblastx{
      input:
        Fname      = fname,
        Queryfa    = queryfa,
        docker     = blast_docker,
        outfmt     = outfmt,
        evalue     = evalue,
        Outfile    = Outfile,
        threads    = threads
    }
  }
  output {
    File fina_output =select_first([runtblastx.out,runblastp.out,runblastn.out,runblastx.out,runtblastn.out])
  }
}

task runblastn {
  input {
    String  docker
    File    Queryfa
    String  Fname
    String  Outfile
    Int threads
    #blast optional
    Int     outfmt
    Int     max_target_seqs = 100
    Float   evalue
    Int     word_size       = 28
    Int     reward          = 1
    Int     penalty         = -2
    String  strand          = 'both'
    Int     gapopen         = 0
    Int     gapextend       = 0
    String  dust            = "'20 64 1'"
    Int?    max_hsps
    String  tasks           = "megablast"
    String?   taxids
    String? negative_taxids
    Boolean lcase_masking   = false
  }
  runtime{
    docker : docker
    cpu    : "8"
    memory : "16G"
  }
  command {
    set -e
    blastn -db "${Fname}" \
    -show_gis \
    -query ${Queryfa} \
    -outfmt ${outfmt} \
    -out     ${Outfile} \
    -max_target_seqs ${max_target_seqs} \
    -evalue ${evalue} \
    -word_size ${word_size} \
    -penalty ${penalty} \
    -reward  ${reward} \
    -dust ${dust} \
    -gapopen ${gapopen} \
    -gapextend ${gapextend} \
    -task ${tasks} \
    -strand  ${strand} \
    -num_threads ${threads} \
    ${true='-lcase_masking' false='' lcase_masking} ${"-max_hsps "+max_hsps} ${"-taxids " +taxids} ${"-negative_taxids " +negative_taxids}\

  }
  output {
    File out = "${Outfile}"
  }
}

task runblastp {
  input {
    String  docker
    File    Queryfa
    String  Fname
    #blast optional
    Int     outfmt
    String  Outfile
    Float   evalue
    Int threads
    Int     max_target_seqs  = 100
    Int     word_size        = 6
    String seg               = "no"
    String comp_based_stats  = "2"
    String matrix            = "BLOSUM62"
    Int     gapopen          = 11
    Int     gapextend        = 1
    Int?    max_hsps
    String?   taxids
    String? negative_taxids
    Boolean lcase_masking    = false
  }
  runtime{
    docker : docker
    cpu    : "8"
    memory : "16G"
  }
  command {
    set -e
    blastp -db "${Fname}" \
    -query ${Queryfa} \
    -outfmt ${outfmt} \
    -out     ${Outfile} \
    -max_target_seqs ${max_target_seqs} \
    -comp_based_stats ${comp_based_stats} \
    -evalue ${evalue} \
    -word_size ${word_size} \
    -matrix   ${matrix} \
    -seg     ${seg} \
    -gapopen ${gapopen} \
    -gapextend ${gapextend} \
    -num_threads ${threads} \
    ${true='-lcase_masking' false='' lcase_masking} ${"-max_hsps "+max_hsps} ${"-taxids " +taxids} ${"-negative_taxids " +negative_taxids} \

  }
  output {
    File out = "${Outfile}"
  }
}

task runblastx {
  input {
    File   Queryfa
    String Fname
    Int    outfmt
    Float  evalue
    String Outfile
    String docker
    Int     threads
    Int    max_target_seqs  =  100
    Int    word_size        =  6
    String seg              =  "'12 2.2 2.5'"
    String comp_based_stats =  "2"
    String matrix           =  "BLOSUM62"
    Int    gapopen          =  11
    Int    gapextend        =  1
    String?   taxids
    String? negative_taxids
    Int?    max_hsps
    Boolean lcase_masking    = false
  }
  runtime{
    docker : docker
    cpu   : "8"
    memory : "16G"
  }
  command {
    set -e
    blastx -db "${Fname}" \
    -query ${Queryfa} \
    -outfmt ${outfmt} \
    -out     ${Outfile} \
    -max_target_seqs ${max_target_seqs} \
    -comp_based_stats ${comp_based_stats} \
    -evalue ${evalue} \
    -word_size ${word_size} \
    -matrix   ${matrix} \
    -seg     ${seg} \
    -gapopen ${gapopen} \
    -gapextend ${gapextend} \
    -num_threads ${threads} \
    ${true='-lcase_masking' false='' lcase_masking} ${"-max_hsps "+max_hsps} ${"-taxids " +taxids} ${"-negative_taxids " +negative_taxids}\

  }
  output {
    File out = "$${Outfile}"
  }
}

task runtblastn {
  input {
    File   Queryfa
    String Fname
    Int    outfmt
    Float  evalue
    String Outfile
    String docker
    Int    threads
    Int    max_target_seqs  = 100
    Int    word_size        = 6
    String seg              = "'12 2.2 2.5'"
    String comp_based_stats = "2"
    String matrix           = "BLOSUM62"
    Int    gapopen          = 11
    Int    gapextend        = 1
    Boolean lcase_masking   = false
    Int?    max_hsps
    String?   taxids
    String? negative_taxids
  }
  runtime{
    docker :docker
    cpu    : "8"
    memory : "16G"
  }
  command {
    set -e
    tblastn -db "${Fname}" \
    -query ${Queryfa} \
    -outfmt ${outfmt} \
    -out    ${Outfile} \
    -max_target_seqs ${max_target_seqs} \
    -comp_based_stats ${comp_based_stats} \
    -evalue ${evalue} \
    -word_size ${word_size} \
    -matrix   ${matrix} \
    -seg     ${seg} \
    -gapopen ${gapopen} \
    -gapextend ${gapextend} \
    -num_threads ${threads} \
    ${true='-lcase_masking' false='' lcase_masking} ${"-max_hsps "+max_hsps} ${"-taxids " +taxids} ${"-negative_taxids " +negative_taxids}\

  }
  output {
    File out = "${Outfile}"
  }
}

task runtblastx {
  input {
    File   Queryfa
    String Fname
    Int    outfmt
    String Outfile
    Int threads
    Float  evalue
    String docker
    String?   taxids
    Int    word_size        = 3
    Int    max_target_seqs  = 100
    String seg              = "'12 2.2 2.5'"
    String matrix           = "BLOSUM62"
    Boolean lcase_masking   = false
    String? negative_taxids
    Int?    max_hsps
  }
  runtime{
    docker :docker
    cpu    : "8"
    memory : "16G"
  }
  command {
    set -e
    tblastx -db "${Fname}" \
    -query ${Queryfa} \
    -outfmt ${outfmt} \
    -out     ${Outfile} \
    -max_target_seqs ${max_target_seqs} \
    -evalue ${evalue} \
    -word_size ${word_size} \
    -matrix   ${matrix} \
    -seg     ${seg} \
    -num_threads ${threads} \
    ${true='-lcase_masking' false='' lcase_masking} ${"-max_hsps "+max_hsps} ${"-taxids " +taxids} ${"-negative_taxids " +negative_taxids}\
  }
  output {
    File out = "${Outfile}"
  }
}


