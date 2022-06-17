version 1.0

## Copyright CHUV, 2020
## Script to annotate vcf files with annovar 
## per scattered contigs after jointgenotyping

## Local import
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/JointCalling-tasks-test.wdl" as Tasks

#################################################################
# WORKFLOW DEFINITION - VCF annotation
#################################################################

workflow JointAnnotation {

String pipeline_version = "1.2"

input {
File GATK
File tabix
File conf_file
Array[File] annovar_vcf
String callset_name
}

scatter (idx in range(length(annovar_vcf))) {

    call Tasks.vcfannoScatteredVCF as vcfannoScatteredVCF { 
    input:
    input_vcf = annovar_vcf[idx],
    conf_file = conf_file,
    base_vcf_name = callset_name + "." + idx
    }
}

call Tasks.GatherVcfs as GatherVcfs {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = vcfannoScatteredVCF.output_vcf,
      output_vcf_name = callset_name + ".VEP.annotated.vcf.gz"
  }


  output {
    File final_vcf = GatherVcfs.output_vcf
    File final_vcf_index = GatherVcfs.output_vcf_index
  }

}