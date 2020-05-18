version 1.0


#################################################################
# WORKFLOW DEFINITION
#################################################################

workflow GenotypeRefinement {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File vcf
    File vcf_index
    File one_thousand_genomes_vcf
    File one_thousand_genomes_vcf_index
    String callset_name
  }

call CalculatGenotypePosteriors {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      vcf = vcf,
      vcf_index = vcf_index,
      one_thousand_genomes_vcf = one_thousand_genomes_vcf,
      one_thousand_genomes_vcf_index = one_thousand_genomes_vcf_index,
      output_vcf_filename = callset_name + "refined_priors.vcf.gz"
  }

call VariantFiltration {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      vcf = CalculatGenotypePosteriors.output_vcf,
      vcf_index = CalculatGenotypePosteriors.output_vcf_index,
      output_vcf_filename = callset_name + "refined_priors-filtered.vcf.gz"
  }

  output {
 # multisample VCF file with refined genotypes
  File refined_vcf = CalculatGenotypePosteriors.output_vcf
  File refined_vcf_index = CalculatGenotypePosteriors.output_vcf_index

 # multisample VCF file with refined AND filtered genotypes
  File final_vcf = VariantFiltration.output_vcf
  File final_vcf_index = VariantFiltration.output_vcf_index
  }  
}


#################################################################
# TASKS DEFINITION
#################################################################

# refines genotype quality using all cohort samples (if >10) and the 1000G dataset
task CalculatGenotypePosteriors {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File vcf
    File vcf_index
    File one_thousand_genomes_vcf
    File one_thousand_genomes_vcf_index
    String output_vcf_filename
  }
  command {
    java -Xmx4g \
    -jar ~{GATK} \
    CalculateGenotypePosteriors \
    -R ~{ref_fasta} \
    --supporting ~{one_thousand_genomes_vcf} \
    -V ~{vcf} \
    -O ~{output_vcf_filename}
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "6000"
	queue: "normal"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

# annotates variants with refined GQ <20 as lowQ
task VariantFiltration {
 input {
    File GATK
    File ref_fasta
    File ref_index
    File vcf
    File vcf_index
    String output_vcf_filename
  }
  command {
    java -Xms3g \
    -jar ~{GATK} \
    VariantFiltration \
    -R ~{ref_fasta} \
    -V ~{vcf} \
    -G-filter "GQ < 20.0" \
    -G-filter-name lowGQ \
    -O ~{output_vcf_filename}
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "4000"
	queue: "normal"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}