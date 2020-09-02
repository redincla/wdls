##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-sv/Delly/1/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

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

import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/delly-tasks.wdl" as Tasks

workflow Delly {
  # Run Delly SV detection algorithm on whole genomes in array of bam or cram files.
  input {
#    File full_map # col 1: sample ID , col 2: input bam/cram , col3: input bam/cram index
    Array[File] input_bams #.bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk.
    Array[File] input_bams_index #"[optional] Index for bam_or_cram_file. If not specified, index is assumed to be at bam_file_path + '.bai' or cram_file_path + '.crai'
    Array[String] sample_ids #sample name. Outputs will be sample_name + 'delly.vcf.gz' and sample_name + 'delly.vcf.gz.tbi'
    File ref_fasta #.fasta file with reference used to align bam or cram file
    File ref_index #[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file
    File exclude_regions_bed #text file with lines specifying genomic intervals where SVs should not be called.Each line in (tab-separated) format
    Array[String]? sv_types #[optional] array of event types for Delly to search for. Defaults to ['DEL', 'DUP', 'INV'].
    String cohort_name
  }

    Array[String] event_types = select_first([sv_types, ["DEL", "DUP", "INV"]])
#    Array[Array[String]] inputSamples = read_tsv(full_map)


  scatter (event_type in event_types) {
      scatter (idx in range(length(input_bams))) {  
            call Tasks.gSVCalling as gSVCalling {
                input:
                    input_bam = input_bams[idx],
                    input_bam_index = input_bams_index[idx],
                    sample_id = sample_ids[idx],
                    ref_fasta = ref_fasta,
                    ref_index = ref_index,
                    exclude_regions_bed = exclude_regions_bed,
                    event_type = event_type
            }
      }
    
      call Tasks.MergeBCF as MergeBCF{
      input:
        input_bcfs = gSVCalling.output_bcf,
        cohort_name = cohort_name,
        event_type = event_type
      }

      scatter (idx in range(length(input_bams))) {  
            call Tasks.GenotypeBCF as GenotypeBCF {
                input:
                    input_bam = input_bams[idx],
                    input_bam_index = input_bams_index[idx],
                    sample_id = sample_ids[idx],
                    ref_fasta = ref_fasta,
                    ref_index = ref_index,
                    exclude_regions_bed = exclude_regions_bed,
                    merged_bcf = MergeBCF.merged_bcf,
                    event_type = event_type
            }
      }

      call Tasks.MergeGenotypedBCF as MergeGenotypedBCF {
      input:
        input_bcfs = GenotypeBCF.genotyped_bcf,
        cohort_name = cohort_name,
        event_type = event_type
      }

      call Tasks.FilterGenotypedBCF as FilterGenotypedBCF {
      input:
        input_bcfs = GenotypeBCF.genotyped_bcf,
        cohort_name = cohort_name,
        event_type = event_type
      }
  }

  call Tasks.FinalGatherBCF as FinalGatherBCF {
    input:
      input_bcfs = FilterGenotypedBCF.filtered_genotyped_bcf,
      cohort_name = cohort_name
  }

  output {
    File vcf = FinalGatherBCF.vcf
    File index = FinalGatherBCF.index
  }
}