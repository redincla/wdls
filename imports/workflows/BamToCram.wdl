version 1.0

## Local import
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/Qc.wdl" as QC
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/BamProcessing.wdl" as Processing


#################################################################
# WORKFLOW DEFINITION
#################################################################
workflow BamToCram {
  input {
    File SAMTOOLS
    File PICARD
    File input_bam
    File ref_fasta
    File ref_index
    File ref_dict
    File duplication_metrics
    File chimerism_metrics
    String base_file_name
  }
  # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
  # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
  Float max_duplication_in_reasonable_sample = 0.30
  Float max_chimerism_in_reasonable_sample = 0.15

 ### [1] Convert the final merged recalibrated BAM file to CRAM format
  call Processing.ConvertToCram as ConvertToCram {
    input:
        SAMTOOLS = SAMTOOLS,
        input_bam = input_bam,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        output_basename = base_file_name
  }

 ### [2] Check whether the data has massively high duplication or chimerism rates
  call QC.CheckPreValidation as CheckPreValidation {
    input:
      duplication_metrics = duplication_metrics,
      chimerism_metrics = chimerism_metrics,
      max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample
 }

### [3] Validate the CRAM file
  call QC.ValidateSamFile as ValidateCram {
    input:
        PICARD = PICARD,
        input_bam = ConvertToCram.output_cram,
        input_bam_index = ConvertToCram.output_cram_index,
        report_filename = base_file_name + ".cram.validation_report",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ignore = ["MISSING_TAG_NM"],
        max_output = 1000000000,
        is_outlier_data = CheckPreValidation.is_outlier_data
  }

  output {
     File output_cram = ConvertToCram.output_cram
     File output_cram_index = ConvertToCram.output_cram_index
     File output_cram_md5 = ConvertToCram.output_cram_md5
     File validate_cram_file_report = ValidateCram.report
  }
}