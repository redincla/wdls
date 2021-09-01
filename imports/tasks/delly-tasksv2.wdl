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
  }

  Boolean is_bam = basename(input_bam, ".bam") + ".bam" == basename(input_bam)
  # ensure there's sufficient memory. These coefficients are obtained
  # via linear regression to test data, assuming uncertainty in memory
  # usage proportional to bam or cram size, and adding 3 standard
  # deviations to best fit estimate of memory usage.
  Float bam_or_cram_size = size(input_bam, "GiB")
  Boolean bam_is_big = if (bam_or_cram_size > 10) then true else false
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

#-q 20 -s 15 options for large bam files, else stuck. Also possible, tweek r: -r 50 or to exclude all translocation calls -r 200
  command <<<
    module add UHTS/Analysis/delly/0.7.8
    export OMP_NUM_THREADS=2
    BCF="~{sample_id}.delly.bcf"
    delly call \
      -g "~{ref_fasta}" \
      -x "~{exclude_regions_bed}" \
      -o "$BCF" \
      -n "~{input_bam}" \
      ~{true='-q 20 -s 15 -r 50' false='' bam_is_big} \ 
  >>>

  output {
    File output_bcf = "~{sample_id}.delly.bcf"
  }

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "15000"  #10000 for large WGS initially ~{mem_size_Mb}, trying to increase for stuck samples
  runtime_minutes: "10000"
  }
}

############
### Merge bcf 
############

task MergeBCF {
  input {
    Array[File] input_bcfs
    String base_name
  }
  command <<<
    module add UHTS/Analysis/delly/0.7.8
    export OMP_NUM_THREADS=2
    BCF="~{base_name}.delly.bcf"
    delly merge ~{sep=' ' input_bcfs} -o "$BCF"
  >>>

  output {
    File merged_bcf = "~{base_name}.delly.bcf"
  }

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "8000"
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
  }
  command <<<
    module add UHTS/Analysis/delly/0.7.8
    export OMP_NUM_THREADS=2
    BCF="~{sample_id}.delly.geno.bcf"
    delly call -g ~{ref_fasta} -v ~{merged_bcf} -o "$BCF" -x ~{exclude_regions_bed} "~{input_bam}"
  >>>

  output {
    File genotyped_bcf = "~{sample_id}.delly.geno.bcf"
    File genotyped_bcf_index = "~{sample_id}.delly.geno.bcf.csi"
  }

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "10000"
  runtime_minutes: "2000"
  }
}

############
### Merge genotyped BCF
############

task MergeGenotypedBCF {
  input {
    Array[File] input_bcfs
    Array[File] input_bcfs_index
    String cohort_name
  }
  command <<<
    module add UHTS/Analysis/samtools/1.10
    BCF="~{cohort_name}.delly.geno.bcf"
    BCF_index="~{cohort_name}.delly.geno.bcf.csi"
    bcftools merge -m id -O b -o "$BCF" ~{sep=' ' input_bcfs} 
    bcftools index "$BCF" -o "$BCF_index"
  >>>

  output {
    File merged_genotyped_bcf = "~{cohort_name}.delly.geno.bcf"
    File merged_genotyped_bcf_index = "~{cohort_name}.delly.geno.bcf.csi"
  }

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "8000"
# runtime_minutes: "580"
  }
}

############
### filter gSVs: at least 20 unrelated samples
############

task FilterGenotypedBCF {
  input {
    File input_bcf
    File input_bcf_index
    String cohort_name
  }
  command <<<
    module add UHTS/Analysis/delly/0.7.8
    module add UHTS/Analysis/EPACTS/3.2.6
    module add UHTS/Analysis/samtools/1.10
    export OMP_NUM_THREADS=2

    BCF="~{cohort_name}.delly.filtered.geno.bcf"
    delly filter -f germline ~{input_bcf} -o "$BCF" 

    VCF_OUT="~{cohort_name}.delly.vcf.gz"
    bcftools convert -O z -o "$VCF_OUT" "$BCF" 
    tabix "$VCF_OUT"
  >>>

  output {
    File final_vcf = "~{cohort_name}.delly.vcf.gz"
    File final_vcf_index = "~{cohort_name}.delly.vcf.gz.tbi"
  }

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "10000"
#    runtime_minutes: "480"
  }
}

############
### Convert BCF to vcf
############

task BCFsToVCFs {
  input {
    File input_bcf
    String base_name
  }

  command <<<
    module add UHTS/Analysis/samtools/1.10
    module add UHTS/Analysis/EPACTS/3.2.6
    module add UHTS/Analysis/vcftools/0.1.15
    set -Eeuo pipefail

    BCF="~{input_bcf}"
    echo "Extracting vcf header"
    bcftools view ~{input_bcf} | grep "^#" > header.vcf
    VCF=$(echo "$BCF" | sed -e 's/.bcf$/.vcf/')
    echo "Converting $BCF to $VCF, skipping header"
    bcftools view $BCF | grep -v "^#" > $VCF
    VCFS="header.vcf $VCF"

    VCF_OUT="~{base_name}.delly.vcf.gz"
    echo "Concatenating vcf into $VCF_OUT"
    cat $VCFS \
      | vcf-sort -c \
      | bcftools reheader -s <(echo "~{base_name}") \
      | bgzip -c \
      > $VCF_OUT
    echo "Indexing $VCF_OUT"
    tabix "$VCF_OUT"
  >>>

  output {
    File vcf = "${base_name}.delly.vcf.gz"
    File index = "${base_name}.delly.vcf.gz.tbi"
  }

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "4000"
#    runtime_minutes: "480"
  }

}