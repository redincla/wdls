version 1.0

## Copyright CHUV, 2021
## Script to streamline JointGenotyping, variant annotation and variant filtering 
## Each workflow could be launched and tuned independently, if needed

## Local import
import "/home/credin/scratch/WGS/wdls/imports/workflows/JointGenotyping.wdl" as JointGenotyping
import "/home/credin/scratch/WGS/wdls/imports/workflows/annotateVCF.wdl" as Annotation
import "/home/credin/scratch/WGS/wdls/imports/workflows/FilterVCF.wdl" as FilterAndSplit

#################################################################
# WORKFLOW DEFINITION - JointCalling to annotated gVCFs
#################################################################

workflow JointGenotypingToAnnotatedVCF {
    String pipeline_version = "1.0"

    input {
    File full_map  #col1: sample name, col2: path to gvcf, col3: path to gvcf index
    File sample_name_map #col1: sample name, col2: path to gvcf. Automatically generated
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
    File one_thousand_genomes_resource_NO_MULTIALLELIC_vcf
    File one_thousand_genomes_resource_NO_MULTIALLELIC_vcf_index

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


    Array[String] sample_list ## deduce from full_map?
    File input_vcf
    File input_vcf_index

    Int cohort_AC_threshold
    Float pop_AF_threshold

    }

    call JointGenotyping.JointGenotyping {
    input:
        GATK = GATK,
        PICARD = PICARD,
        full_map = full_map, 
        base_file_name = base_file_name,
        make_fofn = make_fofn
     }

    call Annotation.JointAnnotation {
    input:
        input_vcf = JointGenotyping.output_vcfs,
        input_vcf_index = JointGenotyping.output_vcf_indices,
        AnnovarDB = AnnovarDB, ##path to AnnovarDB (hg19/hg38)
        genome_build = genome_build, ## hg38 / hg18 / hg19. Needs corresponding dbs to be downloaded!
        GATK = GATK,
        tabix = tabix,
        conf_file = conf_file,
        scatter_count = scatter_count,
        unpadded_intervals_file = unpadded_intervals_file,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict
     }

    call FilterAndSplit.VCFSplitAndFilter {
    input:
        sample_list = sample_list,
        input_vcf = JointAnnotation.final_vcf,
        input_vcf_index = JointAnnotation.final_vcf_index,
        GATK = GATK,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        cohort_AC_threshold = cohort_AC_threshold,
        pop_AF_threshold = pop_AF_threshold
     }

    output {

        File jointgenotyping_detail_metrics_file = JointGenotyping.detail_metrics_file
        File jointgenotyping_metrics_file = JointGenotyping.summary_metrics_file
        File cohort_vcf = JointGenotyping.output_vcfs
        File cohort_vcf_index = JointGenotyping.output_vcf_indices
        File? fingerprint_check = JointGenotyping.crosscheck_fingerprint_check
        Array[File] output_intervals = JointGenotyping.output_intervals

        File cohort_annotated_vcf = JointAnnotation.final_vcf
        File cohort_annotated_vcf_index = JointAnnotation.final_vcf_index

        Array[File] individual_annotated_filtered_vcfs = VCFSplitAndFilter.scattered_vcfs
        Array[File] individual_annotated_filtered_vcfs_index = VCFSplitAndFilter.scattered_vcfs_index
        Array[File] individual_annotated_filtered_tsvs = VCFSplitAndFilter.scattered_tables
          }
}