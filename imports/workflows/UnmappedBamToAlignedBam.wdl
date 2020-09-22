version 1.0

## Copyright Broad Institute, 2019
##
## This WDL defines tasks used for alignment of human whole-genome or exome sequencing data.
##
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

## Local import
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/Alignment.wdl" as Alignment
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/workflows/SplitLargeReadGroup.wdl" as SplitRG
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/Qc.wdl" as QC
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/BamProcessing.wdl" as Processing
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/Utilities.wdl" as Utils


############
### [0.0] Add: FastqToSam https://github.com/jharri34/krueger/issues/9
### or : IlluminaBasecallsToSam from picard: https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/illumina/IlluminaBasecallsToSam.java

#################################################################
# WORKFLOW DEFINITION
#################################################################

workflow UnmappedBamToAlignedBam {
  input {
    Array[File] flowcell_unmapped_bams
    String base_file_name
    File ref_fasta
    File ref_index
    File ref_dict
    File ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File PICARD
    File BWA
    File GATK
    File SAMTOOLS
    File VerifyBamID
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File dbsnp_vcf
    File dbsnp_vcf_index
    Boolean disable_sanity_check

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String recalibrated_bam_basename
  }
    Float cutoff_for_large_rg_in_gb = 20.0
    Int compression_level = 2

    ### [0.1] Retrieve bwa version to include in the PG record in the header of the BAM produced by MergeBamAlignment
    call Alignment.GetBwaVersion as GetBwaVersion { input: BWA = BWA }

    scatter (unmapped_bam in flowcell_unmapped_bams) {
        Float unmapped_bam_size = size(unmapped_bam, "GiB")
        String unmapped_bam_basename = basename(unmapped_bam, ".bam")

        ### [0.2] QC the unmapped BAM with Piccard
        call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics {
        input:
            PICARD = PICARD,
            input_bam = unmapped_bam,
            metrics_filename = unmapped_bam_basename + ".unmapped.quality_yield_metrics"
        }

        ### [0.3] QC the unmapped BAM with fastqc
        call QC.CollectFastQCMetrics as CollectFastQCMetrics {
        input:
            input_bam = unmapped_bam,
            metrics_basename = unmapped_bam_basename
        }
        
        ### [1.0] Align flowcell-level unmapped input bams in parallel
        if (unmapped_bam_size <= cutoff_for_large_rg_in_gb) {
            call Alignment.SamToFastqAndBwaMem as SamToFastqAndBwaMem {
            input:
                input_bam = unmapped_bam,
                PICARD = PICARD,
                BWA = BWA,
                output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
                bwa_version = GetBwaVersion.bwa_version,
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
        }

        ### [1.1] Split bam into multiple smaller bams, map reads to reference and recombine into one bam
        if (unmapped_bam_size > cutoff_for_large_rg_in_gb) {
            call SplitRG.SplitLargeReadGroup as SplitRG {
            input:
                input_bam = unmapped_bam,
                BWA = BWA,
                PICARD = PICARD,
                SAMTOOLS = SAMTOOLS,
                bwa_version = GetBwaVersion.bwa_version,
                output_bam_basename = unmapped_bam_basename + ".aligned.unsorted",
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
        }    

        File output_aligned_bam = select_first([SamToFastqAndBwaMem.output_bam, SplitRG.aligned_bam])

        ### [1.1] QC the aligned but unsorted readgroup BAM
        call QC.CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics {
        input:
            input_bam = output_aligned_bam,
            PICARD = PICARD,
            output_bam_prefix = unmapped_bam_basename + ".readgroup"
        }
    }   

### [2.0] Aggregate aligned+merged flowcell BAM files and mark duplicates
    call Processing.MarkDuplicates as MarkDuplicates {
        input:
        PICARD = PICARD,
        input_bams = output_aligned_bam,
        output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
        metrics_filename = base_file_name + ".duplicate_metrics",
        compression_level = compression_level
    }

### [3.0] Sort aggregated+deduped BAM file and fix tags
    call Processing.SortSam as SortSampleBam {
        input:
        PICARD = PICARD,
        input_bam = MarkDuplicates.output_bam,
        output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
        compression_level = compression_level
    }

### [4.0] Create list of sequences for scatter-gather parallelization
  call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict
  }


### [5.0] Estimate level of cross-sample contamination
  call Processing.CheckContamination as CheckContamination {
    input:
      VerifyBamID = VerifyBamID,  
      input_bam = SortSampleBam.output_bam,
      input_bam_index = SortSampleBam.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      output_prefix = base_file_name + ".preBqsr",
      disable_sanity_check = disable_sanity_check,
      contamination_underestimation_factor = 0.75
  }

### [5.0]  Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    call Processing.BaseRecalibrator as BaseRecalibrator {
      input:
        GATK = GATK,
        input_bam = SortSampleBam.output_bam,
        input_bam_index = SortSampleBam.output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_index = ref_index
    }
  }

 ### [6.0] Merge the recalibration reports resulting from by-interval recalibration
  call Processing.GatherBqsrReports as GatherBqsrReports {
    input:
        GATK = GATK,
        input_bqsr_reports = BaseRecalibrator.recalibration_report,
        output_report_filename = base_file_name + ".recal_data.csv"
  }

### [7.0] Apply the recalibration model by interval
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    call Processing.ApplyBQSR as ApplyBQSR {
      input:
        GATK = GATK,
        input_bam = SortSampleBam.output_bam,
        input_bam_index = SortSampleBam.output_bam_index,
        output_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        compression_level = compression_level
    }
  }

### [8.0] Merge the recalibrated BAM files resulting from by-interval recalibration
  call Processing.GatherSortedBamFiles as GatherBamFiles {
    input:
        PICARD = PICARD,
        input_bams = ApplyBQSR.recalibrated_bam,
        output_bam_basename = base_file_name,
        compression_level = compression_level
  }

output {
    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics
    Array[File] fastqc_metrics_summary = CollectFastQCMetrics.summary_metrics
    Array[File] fastqc_metrics_html = CollectFastQCMetrics.html_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination
    
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    
    File output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File output_bam = GatherBamFiles.output_bam
    File output_bam_index = GatherBamFiles.output_bam_index
  }
}