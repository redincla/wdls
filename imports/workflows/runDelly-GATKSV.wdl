##########################################################################################

## Base script:   https://github.com/broadinstitute/gatk-sv/blob/master/wdl/Delly.wdl
## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements SV calling with DELLY2 by Tobias Rausch (https://github.com/dellytools/delly), per sample
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

workflow Delly {
  input {
    File bam_or_cram_file #.bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk.
    File bam_or_cram_index #"[optional] Index for bam_or_cram_file. If not specified, index is assumed to be at bam_file_path + '.bai' or cram_file_path + '.crai'
    String sample_id #sample name. Outputs will be sample_name + 'delly.vcf.gz' and sample_name + 'delly.vcf.gz.tbi'
    File reference_fasta #.fasta file with reference used to align bam or cram file
    File reference_index #[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file
    File exclude_regions_bed #text file with lines specifying genomic intervals where SVs should not be called.Each line in (tab-separated) format
  }

      call Tasks.gSVCalling as gSVCalling {
        input:
            input_bam = bam_or_cram_file,
            input_bam_index = bam_or_cram_index,
            sample_id = sample_id,
            ref_fasta = reference_fasta,
            ref_index = reference_index,
            exclude_regions_bed = exclude_regions_bed
      }
    
      call Tasks.BCFsToVCFs as BCFsToVCFs {
        input:
          input_bcf = gSVCalling.output_bcf,
          base_name = sample_id
      }

  output {
  File vcf = BCFsToVCFs.vcf
  File index = BCFsToVCFs.index
  }
}