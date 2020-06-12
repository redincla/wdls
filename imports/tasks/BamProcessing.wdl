version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for BAM file processing of human whole-genome or exome sequencing data.
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
### Sort BAM file by coordinate order
############
task SortSam {
  input {
    File PICARD
    File input_bam
    String output_bam_basename
    Int compression_level
  }
  # SortSam spills to disk a lot more because we only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar ~{PICARD} \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "5000"
	queue: "normal"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

############
### Mark duplicate reads to avoid counting non-independent observations
############
task MarkDuplicates {
  input {
    File PICARD
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    #Float total_input_size
    Int compression_level
    Int memory_multiplier = 1
  }
    # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
    # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
    # Float md_disk_multiplier = 3
    # Int disk_size = ceil(md_disk_multiplier * total_input_size) + 20

    Int memory_size = ceil(7.5 * memory_multiplier)
    Int java_memory_size = (memory_size - 2)

  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar ~{PICARD} \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "14000"
	queue: "normal"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

############
### Generate Base Quality Score Recalibration (BQSR) model
############
task BaseRecalibrator {
  input {
    File GATK
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_index
  }
  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms5g \
      -jar ~{GATK} \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "6000"
	queue: "normal"
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

############
### Apply Base Quality Score Recalibration (BQSR) model
############
task ApplyBQSR {
  input {
    File GATK
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_index
    Int compression_level
    Int memory_multiplier = 1
  }

  command {
    java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms3000m \
      -jar ~{GATK} \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      --static-quantized-quals 10 \
      --static-quantized-quals 20 \
      --static-quantized-quals 30 \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "14000"
	queue: "normal"
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_checksum = "~{output_bam_basename}.bam.md5"
  }
}

############
### Merge recalibration reports
############

### adding java arguments to GATK: java -Xmx4G -XX:+PrintGCDetails -jar gatk.jar [program arguments] 
### OR: gatk --java-options "-Xmx4G -XX:+PrintGCDetails" [program arguments]

task GatherBqsrReports {
  input {
    File GATK
    Array[File] input_bqsr_reports
    String output_report_filename
  }
  command {
    java -Xms3000m -jar ~{GATK} \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}
    }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "3500"
	queue: "normal"
  }
 output {
    File output_bqsr_report = "~{output_report_filename}"
  }
}


############
### Combine multiple *sorted* BAM files
############ 
task GatherSortedBamFiles {
  input {
    File PICARD
    Array[File] input_bams
    String output_bam_basename
    Int compression_level
  }
  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar ~{PICARD} \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "3000"
	queue: "normal"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

############
### Combine multiple *unsorted* BAM files
############  
task GatherUnsortedBamFiles {
  input {
    File PICARD
    Array[File] input_bams
    String output_bam_basename
    Int compression_level
  }

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar ~{PICARD} \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "3000"
	queue: "normal"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

############
### Convert BAM file to CRAM format
############  
# Note that reading CRAMs directly with Picard is not yet supported
task ConvertToCram {
  input {
    File SAMTOOLS
    File input_bam
    File ref_fasta
    File ref_index
    String output_basename
  }

  command <<<
    set -e
    set -o pipefail

    ~{SAMTOOLS} view -C -T ~{ref_fasta} ~{input_bam} | \
    tee ~{output_basename}.cram | \
    md5sum | awk '{print $1}' > ~{output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ~{ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ~{output_basename}.cram
  >>>

  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "3000"
	queue: "normal"
  }
  output {
    File output_cram = "~{output_basename}.cram"
    File output_cram_index = "~{output_basename}.cram.crai"
    File output_cram_md5 = "~{output_basename}.cram.md5"
  }
}

############
### Estimate contamination
############  
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
# Here, Broad is handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
task CheckContamination {
  input {
    File VerifyBamID
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File ref_fasta
    File ref_index
    String output_prefix
    Float contamination_underestimation_factor
    Boolean disable_sanity_check
  }
  command <<<
    set -e

    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    ~{VerifyBamID} \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    ~{true="--DisableSanityCheck" false="" disable_sanity_check} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('~{output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX(alpha)"])/~{contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
	cpus: "2"
	requested_memory_mb_per_core: "7500"
	queue: "normal"
  }
  output {
    File selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
  }
}

############
### Convert a pair of FASTQs to uBAM
############  
task PairedFastQsToUnmappedBAM {
  input {
    File GATK
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center
  }
  command {
    java -Xmx6g \
    -jar ~{GATK} \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM_UNIT ~{platform_unit} \
    --RUN_DATE ~{run_date} \
    --PLATFORM ~{platform_name} \
    --SEQUENCING_CENTER ~{sequencing_center} 
  }
  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "9000"
    queue: "normal"
  }
  output {
    File output_unmapped_bam = "~{readgroup_name}.unmapped.bam"
  }
}

############
### Validate Bam/Sam file
############  
task ValidateBam {
  input {
    File PICARD
    File input_bam
    File? ref_fasta
    String report_filename 
  }

String ref_param = if defined(ref_fasta) then "REFERENCE_SEQUENCE=~{ref_fasta}" else ""

 command {
   java -Xmx6g -jar ~{PICARD} \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      MODE=SUMMARY 
}
  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "9000"
    queue: "normal"
  }
  output {
    File Bamvalidation_report = "~{report_filename}"
  }
}
