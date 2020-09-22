version 1.0

## Local import
import "/home/credin/scratch/WGS/wdls/imports/tasks/GermlineVariantDiscovery.wdl" as Calling
import "/home/credin/scratch/WGS/wdls/imports/tasks/Qc.wdl" as QC
import "/home/credin/scratch/WGS/wdls/imports/tasks/Utilities.wdl" as Utils
import "/home/credin/scratch/WGS/wdls/imports/tasks/BamProcessing.wdl" as BamProcessing

#################################################################
# WORKFLOW DEFINITION
#################################################################
workflow VariantCalling {
  input {
    File PICARD
    File SAMTOOLS
    File BWA
    File GATK
    File GATK3
    File calling_interval_list
    File evaluation_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Float? contamination
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_index
    File ref_dict
    File dbsnp_vcf
    File dbsnp_vcf_index
    String base_file_name
    String final_vcf_base_name
    Boolean make_gvcf = true
    Boolean make_bamout
    Boolean use_gatk3_haplotype_caller = false
  }
  #parameter_meta {
  #  make_bamout: "For CNNScoreVariants to run with a 2D model, a bamout must be created by HaplotypeCaller. The bamout is a bam containing information on how HaplotypeCaller remapped reads while it was calling variants. See https://gatkforums.broadinstitute.org/gatk/discussion/5484/howto-generate-a-bamout-file-showing-how-haplotypecaller-has-remapped-sequence-reads for more details."
  #}

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utils.ScatterIntervalList as ScatterIntervalList {
    input:
      PICARD = PICARD,
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }

  # Call variants in parallel over WGS calling intervals
  scatter (scattered_interval_list in ScatterIntervalList.out) {

    if (use_gatk3_haplotype_caller) {
      call Calling.HaplotypeCaller_GATK35_GVCF as HaplotypeCallerGATK3 {
        input:
        GATK = GATK,
        GATK3 = GATK3,
        input_bam = input_bam,
        interval_list = scattered_interval_list,
        gvcf_basename = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        contamination = contamination
      }
    }

    if (!use_gatk3_haplotype_caller) {

      # Generate GVCF by interval
      call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
        input:
          GATK = GATK,
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = scattered_interval_list,
          vcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_index = ref_index,
          contamination = contamination,
          make_gvcf = make_gvcf,
          make_bamout = make_bamout
       }


      # If bamout files were created, we need to sort and gather them into one bamout
      if (make_bamout) {
        call BamProcessing.SortSam as SortBamout {
          input:
            PICARD = PICARD,
            input_bam = HaplotypeCallerGATK4.bamout,
            output_bam_basename = final_vcf_base_name,
            compression_level = 2
        }
      }
    }

    File vcfs_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf, HaplotypeCallerGATK4.output_vcf])
    File vcf_indices_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf_index, HaplotypeCallerGATK4.output_vcf_index])
  }

  # Combine by-interval (g)VCFs into a single sample (g)VCF file
  String merge_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  call Calling.MergeVCFs as MergeVCFs {
    input:
      PICARD = PICARD,
      input_vcfs = vcfs_to_merge,
      input_vcfs_indexes = vcf_indices_to_merge,
      output_vcf_name = final_vcf_base_name + merge_suffix
  }

  if (make_bamout) {
    call MergeBamouts {
      input:
        SAMTOOLS = SAMTOOLS,
        bams = select_all(SortBamout.output_bam),
        output_base_name = final_vcf_base_name
    }
  }

  # Validate the (g)VCF output of HaplotypeCaller
  call QC.ValidateVCF as ValidateVCF {
    input:
      GATK = GATK,
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      calling_interval_list = calling_interval_list,
      is_gvcf = make_gvcf
  }

  # QC the (g)VCF
  call QC.CollectVariantCallingMetrics as CollectVariantCallingMetrics {
    input:
      PICARD = PICARD,
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      metrics_basename = final_vcf_base_name,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      ref_dict = ref_dict,
      evaluation_interval_list = evaluation_interval_list,
      is_gvcf = make_gvcf
  }

  output {
    File vcf_summary_metrics = CollectVariantCallingMetrics.summary_metrics
    File vcf_detail_metrics = CollectVariantCallingMetrics.detail_metrics
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
    File? bamout = MergeBamouts.output_bam
    File? bamout_index = MergeBamouts.output_bam_index
  }
}

# This task is here because merging bamout files using Picard produces an error.
task MergeBamouts {
  input {
    File SAMTOOLS
    Array[File] bams
    String output_base_name
  }
  command {
    ~{SAMTOOLS} merge ~{output_base_name}.bam ~{sep=" " bams}
    ~{SAMTOOLS} index ~{output_base_name}.bam
    mv ~{output_base_name}.bam.bai ~{output_base_name}.bai
  }

  output {
    File output_bam = "~{output_base_name}.bam"
    File output_bam_index = "~{output_base_name}.bai"
  }

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "4000"
  }
}