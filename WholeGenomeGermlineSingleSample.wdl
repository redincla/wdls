version 1.0

## Copyright Broad Institute, 2018
## Piou
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
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
import "./imports/workflows/UnmappedBamToAlignedBam.wdl" as ToBam
import "./imports/workflows/AggregatedBamQC.wdl" as AggregatedQC
import "./imports/tasks/Qc.wdl" as QC
import "./imports/workflows/BamToCram.wdl" as ToCram
import "./imports/tasks/BamProcessing.wdl" as Processing
import "./imports/workflows/VariantCalling.wdl" as ToGvcf


#################################################################
# WORKFLOW DEFINITION
#################################################################
workflow WholeGenomeGermlineSingleSample {
    String pipeline_version = "1.3"

  input {
    Array[File] flowcell_unmapped_bams
    String sample_name
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
    File wgs_coverage_interval_list
    File PICARD
    File BWA
    File GATK
    File GATK3
    File SAMTOOLS
    File VerifyBamID
    Boolean provide_bam_output = false  #output aligned bam+bai in addition to cram+crai files
    Boolean use_gatk3_haplotype_caller = false
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File dbsnp_vcf
    File dbsnp_vcf_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    Boolean disable_sanity_check
    File calling_interval_list
    File evaluation_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Boolean make_bamout #output realigned bam+bai files after initial variant calling step
    File final_vcf_base_name
  }

    # Not overridable:
    Int read_length = 250
    Float lod_threshold = -20.0
    String cross_check_fingerprints_by = "READGROUP"
    String recalibrated_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated"


    call ToBam.UnmappedBamToAlignedBam {
    input:
        flowcell_unmapped_bams = flowcell_unmapped_bams,
        base_file_name = base_file_name,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        PICARD = PICARD,
        SAMTOOLS = SAMTOOLS,
        BWA = BWA,
        GATK = GATK,
        VerifyBamID = VerifyBamID,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        known_indels_sites_indices = known_indels_sites_indices,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        disable_sanity_check = disable_sanity_check,
        contamination_sites_ud = contamination_sites_ud,
        contamination_sites_bed = contamination_sites_bed,
        contamination_sites_mu = contamination_sites_mu,
        recalibrated_bam_basename = recalibrated_bam_basename
     }

    call AggregatedQC.AggregatedBamQC {
    input:
        base_recalibrated_bam = UnmappedBamToAlignedBam.output_bam,
        base_recalibrated_bam_index = UnmappedBamToAlignedBam.output_bam_index,
        base_name = base_file_name,
        recalibrated_bam_base_name = recalibrated_bam_basename,
        PICARD = PICARD,
        ref_fasta = ref_fasta,
	      ref_index = ref_index,
        ref_dict = ref_dict
    }

    call ToCram.BamToCram as BamToCram {
    input:
        SAMTOOLS = SAMTOOLS,
        PICARD = PICARD,
        input_bam = UnmappedBamToAlignedBam.output_bam,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        duplication_metrics = UnmappedBamToAlignedBam.duplicate_metrics,
        chimerism_metrics = AggregatedBamQC.agg_alignment_summary_metrics,
        base_file_name = base_file_name
    }

  # QC the sample WGS metrics (stringent thresholds)
  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
        PICARD = PICARD,
        input_bam = UnmappedBamToAlignedBam.output_bam,
        input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
        metrics_filename = base_file_name + ".wgs_metrics",
        wgs_coverage_interval_list = wgs_coverage_interval_list,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        read_length = read_length
  }

  # QC the sample raw WGS metrics (common thresholds)
  call QC.CollectRawWgsMetrics as CollectRawWgsMetrics {
    input:
        PICARD = PICARD,
        input_bam = UnmappedBamToAlignedBam.output_bam,
        input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
        metrics_filename = base_file_name + ".raw_wgs_metrics",
        wgs_coverage_interval_list = wgs_coverage_interval_list,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        read_length = read_length
  }

  # Initial variant calling on single sample
  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      PICARD = PICARD,
      BWA = BWA,
      GATK = GATK,
      GATK3 = GATK3,
      SAMTOOLS = SAMTOOLS,
      calling_interval_list = calling_interval_list,
      evaluation_interval_list = evaluation_interval_list,
      haplotype_scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      contamination = UnmappedBamToAlignedBam.contamination,
      input_bam = UnmappedBamToAlignedBam.output_bam,
      input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      base_file_name = base_file_name,
      final_vcf_base_name = final_vcf_base_name,
      make_bamout = make_bamout
  }

  if (provide_bam_output) {
    File provided_output_bam = UnmappedBamToAlignedBam.output_bam
    File provided_output_bam_index = UnmappedBamToAlignedBam.output_bam_index
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = UnmappedBamToAlignedBam.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = UnmappedBamToAlignedBam.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = UnmappedBamToAlignedBam.unsorted_read_group_quality_distribution_metrics
    
    File selfSM = UnmappedBamToAlignedBam.selfSM
    Float contamination = UnmappedBamToAlignedBam.contamination   

    File read_group_alignment_summary_metrics = AggregatedBamQC.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = AggregatedBamQC.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = AggregatedBamQC.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = AggregatedBamQC.read_group_gc_bias_summary_metrics

    File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = AggregatedBamQC.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = AggregatedBamQC.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = AggregatedBamQC.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = AggregatedBamQC.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = AggregatedBamQC.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
    File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

    File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics

    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    File output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

    File gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    File gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics

    File? output_bam = provided_output_bam
    File? output_bam_index = provided_output_bam_index

    File output_cram = BamToCram.output_cram
    File output_cram_index = BamToCram.output_cram_index
    File output_cram_md5 = BamToCram.output_cram_md5

    File validate_cram_file_report = BamToCram.validate_cram_file_report

    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index
  }
}
