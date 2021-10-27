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
import "/home/credin/scratch/WGS/wdls/imports/tasks/BamProcessing.wdl" as Processing

#################################################################
# WORKFLOW DEFINITION
#################################################################
workflow ConvertPairedFastQsToUnmappedBamWf {
  input {
    File GATK
    String sample_name 
    String fastq_1 
    String fastq_2 
    String readgroup_name 
    String library_name 
    String platform_unit 
    String run_date 
    String platform_name 
    String sequencing_center 
    Boolean make_fofn
  }
  # Convert pair of FASTQs to uBAM
  call Processing.PairedFastQsToUnmappedBAM as PairedFastQsToUnmappedBAM {
    input:
      GATK = GATK,
      sample_name = sample_name,
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,
      readgroup_name = readgroup_name,
      library_name = library_name,
      platform_unit = platform_unit,
      run_date = run_date,
      platform_name = platform_name,
      sequencing_center = sequencing_center,
  }

  #Create a file with the generated ubam
  if (make_fofn) {  
    call CreateFoFN {
      input:
        ubam = PairedFastQsToUnmappedBAM.output_unmapped_bam,
        fofn_name = sample_name + ".ubam"
    }
  }
  
  # Outputs that will be retained when execution is complete
  output {
    File output_unmapped_bam = PairedFastQsToUnmappedBAM.output_unmapped_bam
    File? unmapped_bam_list = CreateFoFN.fofn_list
  }
}

# Creates a file listing paths to all uBAMs (each row = path to a uBAM file)
# In this case there will only be one file path in the txt file but this format is used by 
# the pre-processing for variant discvoery workflow. 
task CreateFoFN {
  input {
    String ubam
    String fofn_name
  }
  command {
    echo ~{ubam} > ~{fofn_name}.list
  }
  output {
    File fofn_list = "~{fofn_name}.list"
  }
  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "1000"
  }
}