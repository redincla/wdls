version 1.0

## workflow to test COllectRawWGS metrics that keep failing through Cromwell


#################################################################
# WORKFLOW DEFINITION
#################################################################

workflow CollectMetrics {
 input {
    File PICARD
    File input_bam
    File input_bam_index
    String base_file_name
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_index
 }

    Int read_length = 250

### [1] Collect Raw Metrics
  call CollectRawWgsMetrics {
    input:
        PICARD = PICARD,
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        metrics_filename = base_file_name + ".wgs_metrics",
        wgs_coverage_interval_list = wgs_coverage_interval_list,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        read_length = read_length
  }

 ### [2] Collect WGS Metrics
 # call CollectWgsMetrics {
 #   input:
 #       PICARD = PICARD,
 #       input_bam = input_bam,
 #       input_bam_index = input_bam_index,
 #       metrics_filename = base_file_name + ".raw_wgs_metrics",
 #       wgs_coverage_interval_list = wgs_coverage_interval_list,
 #       ref_fasta = ref_fasta,
 #       ref_index = ref_index,
 #       read_length = read_length
 #}

 output {
 #   File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics
  }

}

#################################################################
# TASKS DEFINITION
#################################################################


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
    java -Xms2000m -jar ~{PICARD} \
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
	  requested_memory_mb_per_core: "3000"
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
    Int memory_multiplier = 1
  }
    Float ref_size = size(ref_fasta, "GiB") + size(ref_index, "GiB")
    Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

    Int memory_size = ceil((if (disk_size < 110) then 5 else 7) * memory_multiplier)
    String java_memory_size = (memory_size - 1) * 1000

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