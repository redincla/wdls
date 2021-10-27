version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements a split of large readgroups for human whole-genome and exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

## Local import
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/Alignment.wdl" as Alignment
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/BamProcessing.wdl" as Processing

#################################################################
# WORKFLOW DEFINITION
#################################################################
workflow SplitLargeReadGroup {
  input {
        File input_bam
        File BWA
        File PICARD
        File SAMTOOLS
        String bwa_version
        String output_bam_basename
        File ref_fasta
        File ref_index
        File ref_dict
        File ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        Int compression_level
        Int reads_per_file = 48000000
  }
  call Alignment.SamSplitter as SamSplitter {
    input :
      PICARD = PICARD,
      SAMTOOLS = SAMTOOLS,
      input_bam = input_bam,
      n_reads = reads_per_file,
      compression_level = compression_level
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    Float current_unmapped_bam_size = size(unmapped_bam, "GiB")
    String current_name = basename(unmapped_bam, ".bam")

    call Alignment.SamToFastqAndBwaMem as SamToFastqAndBwaMem {
      input:
        input_bam = unmapped_bam,
        PICARD = PICARD,
        BWA = BWA,
        output_bam_basename = current_name + ".aligned.unsorted",
        bwa_version = bwa_version,
        compression_level = compression_level,
        ref_fasta = ref_fasta,
	      ref_index = ref_index,
		    ref_dict = ref_dict,
		    ref_alt = ref_alt,
		    ref_bwt = ref_bwt,
		    ref_amb = ref_amb,
		    ref_ann = ref_ann,
		    ref_pac = ref_pac,
		    ref_sa = ref_sa
    }

    Float current_mapped_size = size(SamToFastqAndBwaMem.output_bam, "GiB")
  }

  call Processing.GatherUnsortedBamFiles as GatherMonolithicBamFile {
    input:
        PICARD = PICARD,
        input_bams = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = output_bam_basename,
        compression_level = compression_level
  }
  output {
    File aligned_bam = GatherMonolithicBamFile.output_bam
  }
}
