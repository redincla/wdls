version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for alignment of human whole-genome or exome sequencing data.
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

# Local Import
# None

############
### Get version of BWA
############
task GetBwaVersion {
  input {
    File BWA
  }
    command {
        ~{BWA} 2>&1 | \
        grep -e '^Version' | \
        sed 's/Version: //'
    }
    runtime {
	    runtime_minutes: "30"
	    cpus: "1"
	    requested_memory_mb_per_core: "1000"
	    queue: "normal"
    }
    output {
    String bwa_version = read_string(stdout())
  }
}

############
### Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
############
task SamToFastqAndBwaMem {
  input {
        File input_bam
        File PICARD
        File BWA
        String bwa_version
        String output_bam_basename
        File ref_fasta
        File ref_index
        File ref_dict
        File ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa

        Int compression_level
  }      

    command <<<
    set -o pipefail
    set -e

    # if reference_fasta.ref_alt has data in it,
    if [ -s ~{ref_alt} ]; then
      java -Xms1000m -Xmx1000m -jar ~{PICARD} \
        SamToFastq \
        INPUT=~{input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      ~{BWA} mem -K 100000000 -p -v 3 -t 16 -Y ~{ref_fasta} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx1000m -jar ~{PICARD} \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        REFERENCE_SEQUENCE=~{ref_fasta} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="~{bwa_version}" \
        PROGRAM_GROUP_COMMAND_LINE="~{BWA} mem -K 100000000 -p -v 3 -t 16 -Y ~{ref_fasta}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

        #grep -m1 "read .* ALT contigs" ~{output_bam_basename}.bwa.stderr.log | grep -v "read 0 ALT contigs"

    # else if reference_fasta.ref_alt is empty or could not be found
    else
      exit 1;
    fi
  >>>
    runtime {
	    runtime_minutes: "3000"
	    cpus: "3"
	    requested_memory_mb_per_core: "14000"
	    queue: "normal"
    }
    output {
    File output_bam = "~{output_bam_basename}.bam"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
  }
}

############
### Split bam file when too big
############
task SamSplitter {
  input {
    File PICARD
    File SAMTOOLS
    File input_bam
    Int n_reads
    Int compression_level
    Int? total_reads
  }
  command {
    set -e
    mkdir output_dir

    total_reads=$(~{SAMTOOLS} view -c ~{input_bam}) 

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar ~{PICARD} SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  runtime {
	  runtime_minutes: "8000"
	  cpus: "1"
	  requested_memory_mb_per_core: "3750"
	  queue: "normal"
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
}
