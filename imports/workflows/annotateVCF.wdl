version 1.0

## Copyright CHUV, 2020
## Script to annotate vcf files with annovar 
## per scattered contigs after jointgenotyping

## Local import
import "/home/credin/scratch/WGS/wdls/imports/tasks/JointCalling-tasks-test.wdl" as Tasks

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
File conf_file
}

String callset_name  = basename(input_vcfs[0], ".0.vcf.gz")

scatter (idx in range(length(input_vcfs))) {  
    call Tasks.AnnovarScatteredVCF as AnnovarScatteredVCF { 
    input:
    input_vcf = input_vcfs[idx],
    AnnovarDB = AnnovarDB,
    genome_build = genome_build,
    base_vcf_name = callset_name + "." + idx
    }

    call Tasks.vcfannoScatteredVCF as vcfannoScatteredVCF { 
    input:
    input_vcf = AnnovarScatteredVCF.output_vcf,
    conf_file = conf_file,
    base_vcf_name = callset_name + "." + idx
    }
}

call Tasks.GatherVcfs as GatherVcfs {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = vcfannoScatteredVCF.output_vcf,
      output_vcf_name = callset_name + ".annotated.vcf.gz"
  }


  output {
    Array[File] scattered_vcfs = vcfannoScatteredVCF.output_vcf
    File final_vcf = GatherVcfs.output_vcf
    File final_vcf_index = GatherVcfs.output_vcf_index
  }

}