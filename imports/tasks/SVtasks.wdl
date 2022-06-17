version 1.0

############
### Fix DT field in BAM header to avoid go/hts issues, as called by smoove
############
task FixHeaderCRAM {
  input {
    File input_cram
    File input_cram_index
  }

String base_name = basename(input_cram, ".cram")

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc/9.3.0
    module load samtools/1.12

    cp ~{input_cram} .  #to avoid messing up with original cram file
    samtools view -H ~{input_cram} >> header
    cat header | grep "DT:" | cut -f 8 | sort -u > dates 
    while read line; do
        year=$(echo $line | cut -f 2 -d: | cut -f 1 -d- | cut -c1-2 )
        month=$(echo $line | cut -f 2 -d: | cut -f 1 -d- | cut -c3-4)
        day=$(echo $line | cut -f 2 -d: | cut -f 1 -d- | cut -c5-6)
        DT_fixed="20${year}-${month}-${day}T010000"
        sed -i "s@${line}@DT:${DT_fixed}@g" header
    done < dates

    samtools reheader -P -i header ~{base_name}.cram
    samtools index ~{base_name}.cram
>>>

  output {
    File output_cram = "~{base_name}.cram"
    File output_cram_index = "~{base_name}.cram.crai"
  }

 runtime {
    cpus: "1"
    requested_memory_mb_per_core: "2000"
    runtime_minutes: "120"
  }
}

############
### Fix DT field in BAM header to avoid go/hts issues, as called by smoove
############
task FixHeaderBAM {
  input {
    File input_bam
    File input_bam_index
  }

String base_name = basename(input_bam, ".bam")

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc/9.3.0
    module load samtools/1.12

    cp ~{input_bam} .  #to avoid messing up with original cram file
    samtools view -H ~{input_bam} >> header
    cat header | grep "DT:" | cut -f 8 | sort -u > dates 
    while read line; do
        year=$(echo $line | cut -f 2 -d: | cut -f 1 -d- | cut -c1-2 )
        month=$(echo $line | cut -f 2 -d: | cut -f 1 -d- | cut -c3-4)
        day=$(echo $line | cut -f 2 -d: | cut -f 1 -d- | cut -c5-6)
        DT_fixed="20${year}-${month}-${day}T010000"
        sed -i "s@${line}@DT:${DT_fixed}@g" header
    done < dates

    samtools reheader -P header ~{base_name}.bam  > tmp.bam #### cannot replace the header in-place for bam files... 
    mv tmp.bam ~{base_name}.bam
    samtools index ~{base_name}.bam
>>>

  output {
    File output_bam = "~{base_name}.bam"
    File output_bam_index = "~{base_name}.bam.bai"
  }

 runtime {
    cpus: "1"
    requested_memory_mb_per_core: "2000"
    runtime_minutes: "120"
  }
}

############
# Collect coverage counts  -> double check the  disabled_read_filters_arr input + -Xmx~{command_mem_mb}m mem value
############
task CollectCounts {
  input {
    File intervals
    File bam
    File bam_idx
    String sample_id
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
  }

  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR="/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar"

    gatk --java-options "-Xmx10g" CollectReadCounts \
      -L ~{intervals} \
      --input ~{bam} \
      --read-index ~{bam_idx} \
      --reference ~{ref_fasta} \
      --format TSV \
      --interval-merging-rule OVERLAPPING_ONLY \
      --output ~{sample_id}.counts.tsv 

    sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:~{sample_id}/g" ~{sample_id}.counts.tsv
    bgzip ~{sample_id}.counts.tsv
  >>>

 runtime {
    cpus: "1"
    requested_memory_mb_per_core: "12000"
    runtime_minutes: "720"
  }

  output {
    File counts = "~{sample_id}.counts.tsv.gz"
  }
}


############
# Task to run collect-pesr on a single sample
############
task PESRCollection {
  input {
    File cram
    File cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    String sample_id
  }

  command <<<

    set -euo pipefail

    export GATK_LOCAL_JAR="/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar"

    gatk --java-options "-Xmx4g" PairedEndAndSplitReadEvidenceCollection \
        -I ~{cram} \
        --pe-file ~{sample_id}.disc.txt.gz \
        --sr-file ~{sample_id}.split.txt.gz \
        --sample-name ~{sample_id} \
        -R ~{reference_fasta}

    tabix -f -s1 -b 2 -e 2 ~{sample_id}.disc.txt.gz
    tabix -f -s1 -b 2 -e 2 ~{sample_id}.split.txt.gz

  >>>

   runtime {
    cpus: "1"
    requested_memory_mb_per_core: "5000"
    runtime_minutes: "720"
  }

    output {
    File split_out = "~{sample_id}.split.txt.gz"
    File split_out_index = "~{sample_id}.split.txt.gz.tbi"
    File disc_out = "~{sample_id}.disc.txt.gz"
    File disc_out_index = "~{sample_id}.disc.txt.gz.tbi"
  }
}