version 1.0

##Copyright Broad Institute, 2018
## https://github.com/gatk-workflows/seq-format-conversion/blob/master/paired-fastq-to-unmapped-bam.wdl
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - The following metada descriptors per sample:
##  - readgroup
##  - sample_name
##  - library_name
##  - platform_unit
##  - run_date
##  - platform_name
##  - sequecing_center
##
## Outputs :
## - Set of unmapped BAMs, one per read group
## - File of a list of the generated unmapped BAMs
##
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
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/BamProcessing.wdl" as Processing
import "/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/imports/tasks/Qc.wdl" as QC

#################################################################
# WORKFLOW DEFINITION
#################################################################
workflow pFQtoUnmappedBam {
  input {
    File GATK
    File PICARD
    File full_map # col 1: fastq_1 , col 2: fastq_2 , col3: RG, col4: lib ID, col 5: PU, col6: run date, col7: platform, col8: seq center
    String base_file_name
    Boolean make_fofn
  }

  Array[Array[String]] inputSamples = read_tsv(full_map)

  scatter (line in inputSamples) {
        File fastq_1 = line[0]
        File fastq_2 = line[1]
        String readgroup_name = line[2]
        String library_name = line[3]
        String platform_unit = line[4]
        String run_date = line[5]
        String platform_name = line[6]
        String sequencing_center = line[7]

### [1] Convert pair of FASTQs to uBAM
    call Processing.PairedFastQsToUnmappedBAM as PairedFastQsToUnmappedBAM {
        input:
            GATK = GATK,
            sample_name = base_file_name,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2,
            readgroup_name = readgroup_name,
            library_name = library_name,
            platform_unit = platform_unit,
            run_date = run_date,
            platform_name = platform_name,
            sequencing_center = sequencing_center
    }

#    call Processing.SortSam as SortSampleBam {   useless as automatically query sorted with previous task
#        input:
#          PICARD = PICARD,
#          input_bam = PairedFastQsToUnmappedBAM.output_unmapped_bam,
#          output_bam_basename = sample_name + readgroup_name + ".unaligned_sorted",
#          compression_level = compression_level
#    }

### [2] Validate the BAM file
    call Processing.ValidateBam as ValidateBam {
      input:
          PICARD = PICARD,
          input_bam = PairedFastQsToUnmappedBAM.output_unmapped_bam,
          report_filename = base_file_name + readgroup_name + ".bam.validation_report"
  }

  }
  #Create a file with the generated ubams
  if (make_fofn) {  
    call CreateFoFN {
      input:
        ubam = PairedFastQsToUnmappedBAM.output_unmapped_bam,
        sample_name = base_file_name
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_unmapped_bams = PairedFastQsToUnmappedBAM.output_unmapped_bam
    File? unmapped_bam_list = CreateFoFN.fofn_list
    Array[File] output_bamvalidation_report = ValidateBam.Bamvalidation_report
  }
}

# Creates a file listing paths to all uBAMs (each row = path to a uBAM file)
task CreateFoFN {
  input {
    Array[String] ubam
    String sample_name
  }
  command {
    echo "~{sep='\n' ubam}" >> ~{sample_name}.list
  }
  output {
    File fofn_list = "~{sample_name}.list"
  }
  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "1000"
  }
}