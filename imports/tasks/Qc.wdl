version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for QC of human whole-genome or exome sequencing data.
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

############
### Collect sequencing yield quality metrics
############ 
task CollectQualityYieldMetrics {
  input {
    File PICARD
    File input_bam
    String metrics_filename
  }
  command {
    java -Xms2000m -jar ~{PICARD} \
      CollectQualityYieldMetrics \
      INPUT=~{input_bam} \
      OQ=true \
      OUTPUT=~{metrics_filename}
  }
  runtime {
	runtime_minutes: "120"
	cpus: "1"
	requested_memory_mb_per_core: "3000"
	queue: "normal"
  }
  output {
    File quality_yield_metrics = "~{metrics_filename}"
  }
}

############
### Collect base quality and insert size metrics
############
task CollectUnsortedReadgroupBamQualityMetrics {
  input {
    File PICARD
    File input_bam
    String output_bam_prefix
  }
  command {
    java -Xms5000m -jar ~{PICARD} \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=QualityScoreDistribution \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=ALL_READS

    touch ~{output_bam_prefix}.insert_size_metrics
    touch ~{output_bam_prefix}.insert_size_histogram.pdf
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "7000"
	queue: "normal"
  }
  output {
    File base_distribution_by_cycle_pdf = "~{output_bam_prefix}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "~{output_bam_prefix}.base_distribution_by_cycle_metrics"
    File insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    File quality_by_cycle_pdf = "~{output_bam_prefix}.quality_by_cycle.pdf"
    File quality_by_cycle_metrics = "~{output_bam_prefix}.quality_by_cycle_metrics"
    File quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
  }
}

############
### Collect alignment summary and GC bias quality metrics
############  
task CollectReadgroupBamQualityMetrics {
  input {
    File PICARD
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_index
    Boolean collect_gc_bias_metrics = true
  }
  command {
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics

    java -Xms5000m -jar ~{PICARD} \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "7000"
	queue: "normal"
  }
  output {
    File alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    File gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
  }
}

############
### Collect quality metrics from the aggregated bam
############   
task CollectAggregationMetrics {
  input {
    File PICARD
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_index
    Boolean collect_gc_bias_metrics = true
  }
  command {
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics \
      ~{output_bam_prefix}.insert_size_metrics \
      ~{output_bam_prefix}.insert_size_histogram.pdf

    java -Xms5000m -jar ~{PICARD} \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=QualityScoreDistribution \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "7000"
	queue: "normal"
  }
  output {
    File alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    File bait_bias_detail_metrics = "~{output_bam_prefix}.bait_bias_detail_metrics"
    File bait_bias_summary_metrics = "~{output_bam_prefix}.bait_bias_summary_metrics"
    File gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
    File insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    File pre_adapter_detail_metrics = "~{output_bam_prefix}.pre_adapter_detail_metrics"
    File pre_adapter_summary_metrics = "~{output_bam_prefix}.pre_adapter_summary_metrics"
    File quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
    File error_summary_metrics = "~{output_bam_prefix}.error_summary_metrics"
  }
}

############
### Check whether the data has massively high duplication or chimerism rates
############   
task CheckPreValidation {
input {
    File duplication_metrics
    File chimerism_metrics
    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
  }
  command <<<
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ~{duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ~{chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

    python <<CODE

    import csv
    with open('duplication.csv') as dupfile:
      reader = csv.DictReader(dupfile, delimiter='\t')
      for row in reader:
        with open("duplication_value.txt","w") as file:
          file.write(row['PERCENT_DUPLICATION'])
          file.close()

    with open('chimerism.csv') as chimfile:
      reader = csv.DictReader(chimfile, delimiter='\t')
      for row in reader:
        with open("chimerism_value.txt","w") as file:
          file.write(row['PCT_CHIMERAS'])
          file.close()

    CODE

  >>>
  runtime {
	  runtime_minutes: "60"
	  cpus: "1"
	  requested_memory_mb_per_core: "2000"
	  queue: "normal"
  }
  output {
    Float duplication_rate = read_float("duplication_value.txt")
    Float chimerism_rate = read_float("chimerism_value.txt")
    Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
  }
}

############
### Validate the CRAM file
############   
task ValidateSamFile {
  input {
    File PICARD
    File input_bam
    File? input_bam_index
    String report_filename
    File ref_dict
    File ref_fasta
    File ref_index
    Int? max_output
    Array[String]? ignore
    Boolean? is_outlier_data
    Int memory_multiplier = 1
  }
    Int memory_size = ceil(7 * memory_multiplier)
    Int java_memory_size = (memory_size - 1) * 1000

  command {
    java -Xms~{java_memory_size}m -jar ~{PICARD} \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      ~{"MAX_OUTPUT=" + max_output} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      ~{default='SKIP_MATE_VALIDATION=false' true='SKIP_MATE_VALIDATION=true' false='SKIP_MATE_VALIDATION=false' is_outlier_data} \
      IS_BISULFITE_SEQUENCED=false
  }
  runtime {
	  cpus: "1"
	  requested_memory_mb_per_core: "10000"
	  queue: "normal"
  }
  output {
    File report = "~{report_filename}"
  }
}

############
### Collect sequencing metrics
############   
# /!\ Will break if read lengths in bam > 250bp
# /!\ Using Fast algorithm will break....

task CollectWgsMetrics {
  input {
    File PICARD
    File input_bam
    File input_bam_index
    String metrics_filename
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_index
    Int read_length
  }
  command {
    java -Xms4000m -jar ~{PICARD} \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=false \
      READ_LENGTH=~{read_length}
  }
  runtime {
	  runtime_minutes: "800"
	  cpus: "1"
	  requested_memory_mb_per_core: "6000"
	  queue: "normal"
  } 
  output {
    File metrics = "~{metrics_filename}"
  }
}

############
### Collect raw WGS metrics (commonly used QC thresholds)
############   
# /!\ Will break if read lengths in bam > 250bp.
# Needs 6G of java memory else will break
# /!\ Using Fast algorithm will break...

task CollectRawWgsMetrics {
  input {
    File PICARD
    File input_bam
    File input_bam_index
    String metrics_filename
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_index
    Int read_length
#    Int memory_multiplier = 1
  }
#    Float ref_size = size(ref_fasta, "GiB") + size(ref_index, "GiB")
#    Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

#    Int memory_size = ceil((if (disk_size < 110) then 5 else 7) * memory_multiplier)
#    String java_memory_size = (memory_size - 1) * 1000

  command {
    java -Xms6000m -jar ~{PICARD} \
      CollectRawWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=false \
      READ_LENGTH=~{read_length}
  }
  runtime {
	  runtime_minutes: "800"
	  cpus: "1"
	  requested_memory_mb_per_core: "14000"
	  queue: "normal"
  } 
  output {
    File metrics = "~{metrics_filename}"
  }
}

############
### Generate a checksum per readgroup
############   
task CalculateReadGroupChecksum {
  input {
    File PICARD
    File input_bam
    File input_bam_index
    String read_group_md5_filename
  }
  command {
    java -Xms1000m -jar ~{PICARD} \
      CalculateReadGroupChecksum \
      INPUT=~{input_bam} \
      OUTPUT=~{read_group_md5_filename}
  }
  runtime {
	  runtime_minutes: "60"
	  cpus: "1"
	  requested_memory_mb_per_core: "2000"
	  queue: "normal"
  }
  output {
    File md5_file = "~{read_group_md5_filename}"
  }
}

############
### Validate a (g)VCF with -gvcf specific validation
############   
task ValidateVCF {
  input {
    File GATK
    File input_vcf
    File input_vcf_index
    File ref_fasta
    File ref_index
    File ref_dict
    File dbsnp_vcf
    File dbsnp_vcf_index
    File calling_interval_list
    Boolean is_gvcf = true
  }
  command {
    java -Xms9000m \
      -jar ~{GATK} \
      ValidateVariants \
      -V ~{input_vcf} \
      -R ~{ref_fasta} \
      -L ~{calling_interval_list} \
      ~{true="-gvcf" false="" is_gvcf} \
      --validation-type-to-exclude ALLELES \
      --dbsnp ~{dbsnp_vcf}
  }
  runtime {
	  cpus: "1"
	  requested_memory_mb_per_core: "12000"
  }
}

############
### Collect variant calling metrics from GVCF output
############  
task CollectVariantCallingMetrics {
  input {
    File PICARD
    File input_vcf
    File input_vcf_index
    String metrics_basename
    File dbsnp_vcf
    File dbsnp_vcf_index
    File ref_dict
    File evaluation_interval_list
    Boolean is_gvcf = true
  }
  command {
    java -Xms2000m -jar ~{PICARD} \
      CollectVariantCallingMetrics \
      INPUT=~{input_vcf} \
      OUTPUT=~{metrics_basename} \
      DBSNP=~{dbsnp_vcf} \
      SEQUENCE_DICTIONARY=~{ref_dict} \
      TARGET_INTERVALS=~{evaluation_interval_list} \
      ~{true="GVCF_INPUT=true" false="" is_gvcf}
  }
  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "3000"
  }
  output {
    File summary_metrics = "~{metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "~{metrics_basename}.variant_calling_detail_metrics"
  }
}

############
### Collect FastQC metrics from FQ/bam files
############  
task CollectFastQCMetrics {
  input {
    File input_bam
    String metrics_basename
    Boolean is_bam = true
  }
  command {
    module add UHTS/Quality_control/fastqc/0.11.7
    fastqc ~{input_bam} ~{true="-f bam" false="-f fastq" is_bam}
  }
  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "8000"
    runtime_minutes: "1400"
  }
  output {
    File summary_metrics = "~{metrics_basename}_fastqc.zip"
    File html_metrics = "~{metrics_basename}_fastqc.html"
  }
}