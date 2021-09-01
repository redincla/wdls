version 1.0

## Copyright CHUV, 2021
## Script to reconstruct haplotypes from srWGS data 

#################################################################
# WORKFLOW DEFINITION - HaplotypeSolver
#################################################################

workflow ResolveHaplotye {

String pipeline_version = "1.0"

input {
Array[String] sample_list
File ref_fasta
File ref_index
File ref_dict
File input_vcf
File input_vcf_index
?String target_chromosome
}

String base_name = basename(input_vcf, ".vcf.gz")
#Array[String] sample_list = read_tsv(sample_IDs)
 
scatter (idx in range(length(sample_list))) {  
    call FreeBayes { 
        input:
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            input_bam_or_cram = input_bam_or_cram,
            input_bam_or_cram_index = input_bam_or_cram_index,
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            sample_ID = sample_list[idx],
            base_output_name = base_name,
            target_chromosome = target_chromosome
        }
    }
}

output {
    Array[File] scattered_vcfs = FreeBayes.output_vcf
    Array[File] scattered_vcfs_index = FreeBayes.output_vcf_index
  }

}

#################################################################
# TASK DEFINITION 
#################################################################

############
### Construct haplotypes using freebayes
############
task FreeBayes {
  input {
    File ref_fasta
    File ref_index
    File ref_dict
    File input_bam_or_cram
    File input_bam_or_cram_index
    File input_vcf
    File input_vcf_index
    String base_output_name
    ?String target_chromosome
  }

String chr_arg = if defined(target_chromosome) then "-r ~{target_chromosome}" else ""

command <<<
    module add UHTS/Analysis/freebayes/1.2.0
    module add UHTS/Analysis/EPACTS/3.2.6
    freebayes -f ~{ref_fasta} ~{input_bam_or_cram} -C 5 -g 500 \
    --haplotype-length 50 ~{chr_arg} --haplotype-basis-alleles ~{input_vcf} > ~{base_output_name}_hap.vcf

    bgzip "~{base_output_name}_hap.vcf"
    tabix "~{base_output_name}_hap.vcf.gz"
>>>  

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}_hap.vcf.gz"
    File output_vcf_index = "~{base_output_name}_hap.vcf.gz.tbi"
  }
}