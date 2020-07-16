version 1.0

## Copyright CHUV, 2020
## Script to annotate vcf files with annovar 
## per scattered contigs after jointgenotyping

## Local import
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/JointCalling-tasks-test.wdl" as Tasks

#################################################################
# WORKFLOW DEFINITION - VCF annotation
#################################################################

workflow JointAnnotation {

String pipeline_version = "1.0"

input {
Array[File] input_vcfs
File AnnovarDB   ##path to AnnovarDB (hg19/hg38)
String genome_build ## hg38 / hg18 / hg19. Needs corresponding dbs to be downloaded!
File GATK
File tabix
File bgzip
}

String callset_name  = basename(input_vcfs[0], ".0.vcf.gz")

scatter (idx in range(length(input_vcfs))) {  
    call Tasks.AnnotateScatteredVCF as AnnotateScatteredVCF { ## outputs a vcf and not vcf.gz file!
    input:
    input_vcf = input_vcfs[idx],
    AnnovarDB = AnnovarDB,
    genome_build = genome_build,
    base_vcf_name = callset_name + "." + idx,
    bgzip = bgzip
    }
}

call Tasks.GatherVcfs as GatherVcfs {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = AnnotateScatteredVCF.output_vcf,
      output_vcf_name = callset_name + ".annotated.vcf.gz"
  }


  output {
    Array[File] scattered_vcfs = AnnotateScatteredVCF.output_vcf
    File final_vcf = GatherVcfs.output_vcf
    File final_vcf_index = GatherVcfs.output_vcf_index
  }

}