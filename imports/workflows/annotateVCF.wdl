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
File input_vcf
File input_vcf_index
File AnnovarDB   ##path to AnnovarDB (hg19/hg38)
String genome_build ## hg38 / hg18 / hg19. Needs corresponding dbs to be downloaded!
File GATK
File tabix
File conf_file
Int scatter_count
File unpadded_intervals_file
File ref_fasta
File ref_index
File ref_dict

}

String callset_name  = basename(input_vcf, ".vcf.gz")

call Tasks.SplitIntervalList {
    input:
      GATK = GATK,
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      sample_names_unique_done = true
  }

Array[File] unpadded_intervals = SplitIntervalList.output_intervals

scatter (idx in range(length(unpadded_intervals))) {

    call Tasks.SplitVCF as SplitVCF { 
    input:
    input_vcf = input_vcf,
    input_vcf_index = input_vcf_index,
    interval = unpadded_intervals[idx],
    base_vcf_name = callset_name + "." + idx
    }

    call Tasks.AnnovarScatteredVCF as AnnovarScatteredVCF { 
    input:
    input_vcf = SplitVCF.output_vcf,
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
    File final_vcf = GatherVcfs.output_vcf
    File final_vcf_index = GatherVcfs.output_vcf_index
  }

}