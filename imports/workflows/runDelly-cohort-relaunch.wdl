##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-sv/Delly/1/wdl
## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements SV calling with DELLY2 by Tobias Rausch (https://github.com/dellytools/delly)
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/delly-tasksv2.wdl" as Tasks

workflow DellyCohort {
  # Run Delly SV detection algorithm on whole genomes in array of bam or cram files.
  input {
#    File full_map # col 1: sample ID , col 2: input bam/cram , col3: input bam/cram index
    Array[File] genotyped_bcf 
    Array[File] genotyped_bcf_index 
   String cohort_name
  }

      call Tasks.MergeGenotypedBCF as MergeGenotypedBCF {
      input:
        input_bcfs = genotyped_bcf ,
        input_bcfs_index = genotyped_bcf_index ,
        cohort_name = cohort_name
      }

      call Tasks.FilterGenotypedBCF as FilterGenotypedBCF {
      input:
        input_bcf = MergeGenotypedBCF.merged_genotyped_bcf,
        input_bcf_index = MergeGenotypedBCF.merged_genotyped_bcf_index,
        cohort_name = cohort_name
      }

  output {
    File bcf = FilterGenotypedBCF.final_vcf
    File index = FilterGenotypedBCF.final_vcf_index
  }
}