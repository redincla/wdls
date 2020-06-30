version 1.0

## Copyright Broad Institute, 2019
## 
## This WDL implements the joint discovery and VQSR filtering portion of the GATK 
## Best Practices for germline SNP and Indel discovery in human 
## whole-genome sequencing (WGS) . The workflow requires a sample map 
## file with 50 or more GVCFs and produces a multisample VCF.
##
## Requirements/expectations :
## - One or more GVCFs produced by HaplotypeCaller in GVCF mode 
## - Bare minimum 50 samples. Gene panels are not supported.
##
## Outputs :
## - A VCF file and its index, filtered using variant quality score recalibration 
##   (VQSR) with genotypes for all samples present in the input VCF. All sites that 
##   are present in the input VCF are retained; filtered sites are annotated as such 
##   in the FILTER field.
##
## Note about VQSR wiring :
## The SNP and INDEL models are built in parallel, but then the corresponding 
## recalibrations are applied in series. Because the INDEL model is generally ready 
## first (because there are fewer indels than SNPs) we set INDEL recalibration to 
## be applied first to the input VCF, while the SNP model is still being built. By 
## the time the SNP model is available, the indel-recalibrated file is available to 
## serve as input to apply the SNP recalibration. If we did it the other way around, 
## we would have to wait until the SNP recal file was available despite the INDEL 
## recal file being there already, then apply SNP recalibration, then apply INDEL 
## recalibration. This would lead to a longer wall clock time for complete workflow 
## execution. Wiring the INDEL recalibration to be applied first solves the problem.

## Users working with large sample sets can invoke the GnarlyGenotyper task in the 
## JointGenotyping.wdl workflow. However, the ReblockGVCF tool must be run for all GVCFs 
## produced by HaplotypeCaller before they can be appropriately processed by GnarlyGenotyper. 
## A workflow that applies the reblocking tool is provided here: 
## https://portal.firecloud.org/?return=terra#methods/methodsDev/ReblockGVCF-gatk4_exomes_goodCompression/4

## Cromwell version support 
## - Successfully tested on v47
## - Does not work on versions < v23 due to output syntax
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
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/JointCalling-tasks.wdl" as Tasks

#################################################################
# WORKFLOW DEFINITION - Joint Genotyping for hg38 Whole Genomes and Exomes
#################################################################
workflow JointGenotyping {

  String pipeline_version = "1.2"

input {
  File sample_name_map
  Int sample_num_threshold = 50 #minimum number of samples to get reliable joint calling results

    File GATK
    File tabix
    File PICARD

    File unpadded_intervals_file

    File ref_fasta
    File ref_index
    File ref_dict

    String workspace_dir_name = "genomicsdb"
    String callset_name

    File dbsnp_vcf
    File dbsnp_vcf_index

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69

    Array[String] indel_recalibration_tranche_values
    Array[String] indel_recalibration_annotation_values
    Array[String] snp_recalibration_tranche_values
    Array[String] snp_recalibration_annotation_values

    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf = dbsnp_vcf
    File dbsnp_resource_vcf_index = dbsnp_vcf_index
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index

    Float snp_filter_level
    Float indel_filter_level

    File eval_interval_list

    File haplotype_database

    Int? top_level_scatter_count
    Boolean? gather_vcfs
    Float unbounded_scatter_count_scale_factor = 0.15
    Boolean use_gnarly_genotyper = false  #for large samplesets i.e. >500
    Boolean use_allele_specific_annotations = true
    Boolean cross_check_fingerprints = true
  }

  Boolean allele_specific_annotations = !use_gnarly_genotyper && use_allele_specific_annotations

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000

  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  # We allow overriding this default behavior for testing / special requests.
  Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])

  Int unbounded_scatter_count = select_first([top_level_scatter_count, round(unbounded_scatter_count_scale_factor * num_gvcfs)])
  Int scatter_count = if unbounded_scatter_count > 2 then unbounded_scatter_count else 2 #I think weird things happen if scatterCount is 1 -- IntervalListTools is noop?

  call Tasks.CheckSamplesUnique {
    input:
      sample_name_map = sample_name_map,
      sample_num_threshold = sample_num_threshold
  }

  call Tasks.SplitIntervalList {
    input:
      GATK = GATK,
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      sample_names_unique_done = CheckSamplesUnique.samples_unique
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call Tasks.ImportGVCFs {
      input:
        GATK = GATK,
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        workspace_dir_name = workspace_dir_name, #initially "genomicsdb" but careful if it exists it will be erased!
        batch_size = 50
    }

      call Tasks.GenotypeGVCFs {
        input:
          GATK = GATK,
          workspace_tar = ImportGVCFs.output_genomicsdb,
          interval = unpadded_intervals[idx],
          output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_index = ref_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index
      }

    File genotyped_vcf = GenotypeGVCFs.output_vcf
    File genotyped_vcf_index = GenotypeGVCFs.output_vcf_index

### ToKeep or not? given small sampleset?
    call Tasks.HardFilterAndMakeSitesOnlyVcf {
      input:
        GATK = GATK,
        vcf = genotyped_vcf,
        vcf_index = genotyped_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz"
    }
  }

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      output_vcf_name = callset_name + ".sites_only.vcf.gz"
  }

  call Tasks.IndelsVariantRecalibrator {
    input:
      GATK = GATK,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      mills_resource_vcf = mills_resource_vcf,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = allele_specific_annotations
  }

    call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
        GATK = GATK,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        hapmap_resource_vcf = hapmap_resource_vcf,
        omni_resource_vcf = omni_resource_vcf,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = allele_specific_annotations
    }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    #for really large callsets we give to friends, just apply filters to the sites-only
    call Tasks.ApplyRecalibration {
      input:
        GATK = GATK,
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = SNPsVariantRecalibratorClassic.recalibration,
        snps_recalibration_index = SNPsVariantRecalibratorClassic.recalibration_index,
        snps_tranches = SNPsVariantRecalibratorClassic.tranches,
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        use_allele_specific_annotations = allele_specific_annotations
    }
  }

  # For small callsets we can gather the VCF shards and then collect metrics on it.
    call Tasks.GatherVcfs as FinalGatherVcf {
      input:
        GATK = GATK,
        tabix = tabix,
        input_vcfs = ApplyRecalibration.recalibrated_vcf,
        output_vcf_name = callset_name + ".vcf.gz"
    }

    call Tasks.CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        GATK = GATK,
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict
    }

    if (cross_check_fingerprints) {

      scatter (line in sample_name_map_lines) {
         File gvcf_paths = line[1]
      }

## modified initial GATK version to:  1- add vcf/gvcf indices in the input to save time
##                                    2- use final vcf as 2nd input instead of intermediate recalibrated vcfs
      call Tasks.CrossCheckFingerprint as CrossCheckFingerprintSolo {
        input:
          PICARD = PICARD,
          gvcf_paths = gvcf_paths,
          vcf_paths = ApplyRecalibration.recalibrated_vcf,
          vcf_index_paths = ApplyRecalibration.recalibrated_vcf_index,
          sample_name_map = sample_name_map,
          haplotype_database = haplotype_database,
          output_base_name = callset_name
      }
    }

  # Get the VCFs from either code path
  Array[File?] output_vcf_files = if defined(FinalGatherVcf.output_vcf) then [FinalGatherVcf.output_vcf] else ApplyRecalibration.recalibrated_vcf
  Array[File?] output_vcf_index_files = if defined(FinalGatherVcf.output_vcf_index) then [FinalGatherVcf.output_vcf_index] else ApplyRecalibration.recalibrated_vcf_index

  output {
    # Metrics from either the small or large callset
    File detail_metrics_file = CollectMetricsOnFullVcf.detail_metrics_file
    File summary_metrics_file = CollectMetricsOnFullVcf.summary_metrics_file

    # Outputs from the small callset path through the wdl.
    Array[File] output_vcfs = select_all(output_vcf_files)
    Array[File] output_vcf_indices = select_all(output_vcf_index_files)

    # Output the interval list generated/used by this run workflow.
    Array[File] output_intervals = SplitIntervalList.output_intervals

    # Output the metrics from crosschecking fingerprints.
    File? crosscheck_fingerprint_check = CrossCheckFingerprintSolo.crosscheck_metrics
  }
}