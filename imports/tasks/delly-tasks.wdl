version 1.0

############
### Launch Delly germline caller
############

task gSVCalling {
  input {
    File ref_fasta
    File ref_index
    File exclude_regions_bed
    File input_bam
    File input_bam_index
    String sample_id
    String event_type
  }

  Boolean is_bam = basename(input_bam, ".bam") + ".bam" == basename(input_bam)

  # ensure there's sufficient memory. These coefficients are obtained
  # via linear regression to test data, assuming uncertainty in memory
  # usage proportional to bam or cram size, and adding 3 standard
  # deviations to best fit estimate of memory usage.
  Float bam_or_cram_size = size(input_bam, "GiB")
  Float mem_per_bam_size = 0.03937
  Float mem_bam_offset = 4.4239
  Float mem_per_cram_size = 0.08579
  Float mem_cram_offset = 4.8633
  Float mem_per_read_pairs = "9.745e-9"
  Float mem_read_pairs_offset = 1.947
  Float mem_size_Gb =
      if is_bam then
          mem_per_bam_size * bam_or_cram_size + mem_bam_offset
      else
          mem_per_cram_size * bam_or_cram_size + mem_cram_offset
   Int mem_size_Mb = ceil(mem_size_Gb * 1000)

  command <<<
    module add UHTS/Analysis/delly/0.7.8
    export OMP_NUM_THREADS=2
    BCF="~{sample_id}.delly.~{event_type}.bcf"
    delly call \
      -t ~{event_type} \
      -g "~{ref_fasta}" \
      -x "~{exclude_regions_bed}" \
      -o "$BCF" \
      -n \
      "~{input_bam}"
  >>>

  output {
    File output_bcf = "~{sample_id}.delly.~{event_type}.bcf"
  }

  runtime {
    cpus: "1"
	requested_memory_mb_per_core: "~{mem_size_Mb}"
  }
}

############
### Merge bcf 
############

task MergeBCF {
  input {
    Array[File] input_bcfs
    String cohort_name
    String event_type
  }
  command <<<
    module add UHTS/Analysis/delly/0.7.8
    export OMP_NUM_THREADS=2
    BCF="~{cohort_name}.delly.~{event_type}.bcf"
    delly merge ~{sep=' ' input_bcfs} -o "$BCF"
  >>>

  output {
    File merged_bcf = "~{cohort_name}.delly.~{event_type}.bcf"
  }

  runtime {
    cpus: "1"
	requested_memory_mb_per_core: "5000"
  }
}

############
### genotype BCF across all samples 
############

task GenotypeBCF {
  input {
    File ref_fasta
    File ref_index
    File merged_bcf
    File exclude_regions_bed
    File input_bam
    File input_bam_index
    String sample_id
    String event_type
  }
  command <<<
    module add UHTS/Analysis/delly/0.7.8
    export OMP_NUM_THREADS=2
    BCF="~{sample_id}.delly.~{event_type}.geno.bcf"
    delly call -g ~{ref_fasta} -v ~{merged_bcf} -o "$BCF" -x ~{exclude_regions_bed} "~{input_bam}"
  >>>

  output {
    File genotyped_bcf = "~{sample_id}.delly.~{event_type}.geno.bcf"
  }

  runtime {
    cpus: "1"
	requested_memory_mb_per_core: "5000"
    runtime_minutes: "480"
  }
}

############
### Merge genotyped BCF
############

task MergeGenotypedBCF {
  input {
    Array[File] input_bcfs
    String cohort_name
    String event_type
  }
  command <<<
    module add UHTS/Analysis/samtools/1.10
    BCF="~{cohort_name}.delly.~{event_type}.geno.bcf"
    bcftools merge -m id -O b -o "$BCF" ~{sep=' ' input_bcfs}
  >>>

  output {
    File merged_genotyped_bcf = "~{cohort_name}.delly.~{event_type}.geno.bcf"
  }

  runtime {
    cpus: "1"
	requested_memory_mb_per_core: "5000"
    runtime_minutes: "580"
  }
}

############
### filter gSVs: at least 20 unrelated samples
############

task FilterGenotypedBCF {
  input {
    Array[File] input_bcfs
    String cohort_name
    String event_type
  }
  command <<<
    module add UHTS/Analysis/delly/0.7.8
    export OMP_NUM_THREADS=2
    BCF="~{cohort_name}.delly.~{event_type}.filtered.geno.bcf"
    delly filter -f germline ~{sep=' ' input_bcfs} -o "$BCF" 
  >>>

  output {
    File filtered_genotyped_bcf = "~{cohort_name}.delly.~{event_type}.filtered.geno.bcf"
  }

  runtime {
    cpus: "1"
	requested_memory_mb_per_core: "5000"
    runtime_minutes: "480"
  }
}

############
### Final gather and gvcf conversion
############

task FinalGatherBCF {
  input {
    Array[File] input_bcfs
    String cohort_name
  }

  File first_bcf = input_bcfs[0]
  Int num_bcfs = length(input_bcfs)

  command <<<
    module add UHTS/Analysis/samtools/1.10
    module add UHTS/Analysis/delly/0.7.8
    module add UHTS/Analysis/EPACTS/3.2.6
    export OMP_NUM_THREADS=2
    BCFS="~{sep='\n' input_bcfs}"
    bcftools view ~{first_bcf} | grep "^#" > header.vcf
    VCFS="header.vcf"

    for ((BCF_IND=1; BCF_IND <= ~{num_bcfs}; BCF_IND++)); do
      BCF=$(echo "$BCFS" | sed "$BCF_IND""q;d")
      VCF=$(echo "$BCF" | sed -e 's/.bcf$/.vcf/')
      bcftools view  $BCF | grep -v "^#" > $VCF
      VCFS="$VCFS $VCF"
    done

    VCF_OUT="~{cohort_name}.delly.vcf.gz"
    cat $VCFS \
      | vcf-sort -c \
      | bcftools reheader -s <(echo "~{cohort_name}") \
      | bgzip -c \
      > $VCF_OUT
    tabix "$VCF_OUT"
  >>>

  output {
    File vcf = "${cohort_name}.delly.vcf.gz"
    File index = "${cohort_name}.delly.vcf.gz.tbi"
  }

  runtime {
    cpus: "1"
	requested_memory_mb_per_core: "4000"
    runtime_minutes: "480"
  }
}