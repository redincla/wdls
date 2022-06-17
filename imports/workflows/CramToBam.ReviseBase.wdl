version 1.0

############################################################################################################
### Use to revise Y, R, W, S, K, M, D, H, V, B, X bases in BAM to N for input CRAM files
############################################################################################################

workflow CramToBamReviseBase {
  input {
    File input_cram
    File input_cram_index
    File ref_fasta
    File ref_fasta_fai
    File contiglist
    String base_name
  }

Array[String] contigs = transpose(read_tsv(contiglist))[0]

  scatter (contig in contigs) {
      call SplitCramPerContig {
        input:
          input_cram = input_cram,
          input_cram_index = input_cram_index,
          contig = contig,
          ref_fasta = ref_fasta,
          ref_fasta_fai = ref_fasta_fai,
          base_name = base_name
      }

    call ReviseBaseInBam{
      input:
        bam_file = SplitCramPerContig.bam_file,
        bam_index = SplitCramPerContig.bam_index,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        base_name = base_name
    }
  }
 
  call ConcatBam {
    input:
      base_name = base_name,
      bam_files = ReviseBaseInBam.revised_bam_file,
      bam_indexes = ReviseBaseInBam.revised_bam_index
  }

  output {
    File bam_file = ConcatBam.bam_file
    File bam_index = ConcatBam.bam_index
  }
}


############################################################################################################
#TASK DEFINITION
############################################################################################################

############################################################################################################
#1 split cram files per contig
############################################################################################################
task SplitCramPerContig {
  input {
    File input_cram
    File input_cram_index
    File ref_fasta
    File ref_fasta_fai
    String contig
    String base_name
  }

  command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc/9.3.0
    module load samtools/1.12

        set -Eeuo pipefail
        samtools view \
                 -b \
                 -h \
                 -@ 4 \
                 -T "~{ref_fasta}" \
                 -o "~{base_name}.~{contig}.bam" \
                 "~{input_cram}" \
                 "~{contig}"

        samtools index -@ 4 "~{base_name}.~{contig}.bam"    
  >>>

  output {
    File bam_file  = "~{base_name}.~{contig}.bam"
    File bam_index = "~{base_name}.~{contig}.bam.bai"
  }

 runtime {
    cpus: "4"
    requested_memory_mb_per_core: "1500"
    runtime_minutes: "720"
  }
}

############################################################################################################
#2  Revise base name from cram input file
############################################################################################################

task ReviseBaseInBam{
  input {
    File bam_file
    File bam_index
    File ref_fasta
    File ref_fasta_fai
    String base_name
  }

  command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc/9.3.0
    module load samtools/1.12

    samtools view -H ~{bam_file} > ~{base_name}.revised.sam

    paste \
    <(samtools view ~{bam_file} | cut -f1-9) \
    <(samtools view ~{bam_file} | cut -f10 | sed -e "s/Y/N/g" | sed -e "s/R/N/g" | sed -e "s/W/N/g" | sed -e "s/S/N/g" | sed -e "s/K/N/g" | sed -e "s/M/N/g" | sed -e "s/D/N/g" | sed -e "s/H/N/g" | sed -e "s/V/N/g" | sed -e "s/B/N/g" | sed -e "s/X/N/g" ) \
    <(samtools view ~{bam_file} | cut -f11-) \
    >> ~{base_name}.revised.sam

    samtools view -Sb ~{base_name}.revised.sam -o ~{base_name}.revised.bam

    samtools index ~{base_name}.revised.bam

  >>>

  output{
    File revised_bam_file = "~{base_name}.revised.bam"
    File revised_bam_index = "~{base_name}.revised.bam.bai"
  }

  runtime {
    cpu: "4"
    requested_memory_mb_per_core: "1500"
    runtime_minutes: "720"
  }
}

############################################################################################################
#2  Concat split contig bam files
############################################################################################################
task ConcatBam {
  input {
    String base_name
    Array[File] bam_files
    Array[File] bam_indexes
  }

  command <<<
    samtools merge ~{base_name}.bam ~{sep=" "  bam_files} 
    samtools index ~{base_name}.bam
  >>>

  output{
    File bam_file = "~{base_name}.bam"
    File bam_index = "~{base_name}.bam.bai"
  }

  runtime {
    cpu: "2"
    requested_memory_mb_per_core: "1500"
    runtime_minutes: "720"
  }
}