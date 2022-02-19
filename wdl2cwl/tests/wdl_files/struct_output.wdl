version 1.0

# adapted from https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/pipelines/skylab/build_indices/BuildIndices.wdl

struct References {
  File genome_fa
  File annotation_gtf
}

task GetReferences {
  input {
    String gtf_version
    String organism
    String organism_prefix
  }

  meta {
    description: "Download files needed for building the designated references"
  }

  String ftp_path = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_~{organism}/release_~{gtf_version}"
  String genome_fa = "GRC~{organism_prefix}38.primary_assembly.genome.fa"
  String annotation_gtf = "gencode.v~{gtf_version}.primary_assembly.annotation.gtf"

  command <<<
    set -eo pipefail

    echo a > ~{genome_fa}
    echo b > ~{annotation_gtf}
  >>>

  output {
      References references = object {
      genome_fa: genome_fa,
      annotation_gtf: annotation_gtf
    }
  }

}

task EchoRef {
  input {
    References ref
  }

  String refflat_name = basename(ref.annotation_gtf, ".gtf") + ".refflat.txt"

  command <<<
    touch ~{refflat_name}
  >>>

  output {
   File refflat = refflat_name
  }

}

workflow test {
  input {
    String gtf_version
    String organism
    String organism_prefix
  }

  call GetReferences {
    input: 
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }

  call EchoRef {
    input:
      ref = GetReferences.references
  }

  output {
    File genome_fa = GetReferences.references.genome_fa
    File annotation_gtf = GetReferences.references.annotation_gtf
    File refflat = EchoRef.refflat
  }
  
} 
