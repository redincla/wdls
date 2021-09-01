##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/Manta/1/wdl

##########################################################################################

## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements SV calling with Illumina's Manta software
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

workflow Manta {
  # Run Manta SV detection algorithm on whole genomes from bam or cram files.
  input {
    Array[File] input_bams_or_crams
    Array[File] input_bams_or_crams_indices
    Array[File] sample_ids
    File region_bed
    File region_bed_index
    File samtools
    File GATK
    File ref_fasta #.fasta file with reference used to align bam or cram file
    File ref_index #[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file
    File ref_dict
    String cohort_name
  }

scatter (idx in range(length(input_bams_or_crams))) {  
   call RunManta {
    input:
      bam_or_cram_file = input_bams_or_crams[idx],
      bam_or_cram_index = input_bams_or_crams_indices[idx],
      sample_id = sample_ids[idx],
      reference_fasta = ref_fasta,
      reference_index = ref_index,
      region_bed = region_bed,
      samtools = samtools,
      region_bed_index = region_bed_index
  }
}

call CombineGvCFs {
      input:
        GATK = GATK,
        output_vcf_filename = cohort_name,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        input_gvcf = RunManta.vcf,
        input_gvcf_index = RunManta.index
      }

  output {
    Array[File] vcfs = RunManta.vcf
    Array[File] indices = RunManta.index
  }
}

task RunManta {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File reference_index
    File samtools
    File region_bed
    File region_bed_index
  }

  Boolean is_bam = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  String bam_ext = if is_bam then ".bam" else ".cram"
  String index_ext = if is_bam then ".bai" else ".crai"

  # select number of cpus and jobs
  Int num_cpu_use = 8
  # select number of jobs (threads) to run
  Float jobs_per_cpu_use = 1.3
  Int num_jobs = round(num_cpu_use * jobs_per_cpu_use)

  # ensure there's sufficient memory.
  # NOTE: according to issue #38, manta internally assumes roughly 2GB
  # of RAM per job, although in practice it uses considerably less.
  # 1.5GB is a safe amount
  Float mem_size_gb = num_jobs * 1.5
  Int mem_size_Mb = ceil(mem_size_gb * 1000)
  # ALSO: manta will scale down number of jobs if less than 2GB per
  # job are reported, even if that memory is not needed. The memory
  # reported must be an integer
  Int mem_reported_to_manta = 2 * num_jobs
  
  # ensure there's sufficient disk space
  Float disk_overhead = 10.0
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float bam_or_cram_index_size = size(bam_or_cram_index, "GiB")
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + bam_or_cram_index_size + ref_size + ref_index_size + disk_overhead)

  String expected_index_name = basename(bam_or_cram_file) + index_ext

 command <<<
    module add UHTS/Analysis/EPACTS/3.2.6
    module add UHTS/Analysis/samtools/1.10
    set -Eeuo pipefail

    # if a preemptible instance restarts and runWorkflow.py already
    # exists, manta will throw an error
    if [ -f ./runWorkflow.py ]; then
      rm ./runWorkflow.py
    fi

    ln -s ~{bam_or_cram_file} sample~{bam_ext}
    ln -s ~{bam_or_cram_index} sample~{bam_ext}~{index_ext}


    # prepare the analysis job
    /home/credin/refs/tools/manta-1.6.0.centos6_x86_64/bin/configManta.py \
      --bam sample~{bam_ext} \
      --referenceFasta ~{reference_fasta} \
      --runDir . \
      --callRegions ~{region_bed}

    # always tell manta there are 2 GiB per job, otherwise it will
    # scale back the requested number of jobs, even if they won't
    # need that much memory
    ./runWorkflow.py \
      --mode local \
      --jobs ~{num_jobs} \
      --memGb $((~{num_jobs} * 2))

    # inversion conversion, then compression and index
    python2 /home/credin/refs/tools/manta-1.6.0.centos6_x86_64/libexec/convertInversion.py \
      ~{samtools} \
      ~{reference_fasta} \
      results/variants/diploidSV.vcf.gz \
      | bcftools reheader -s <(echo "~{sample_id}") \
      > diploidSV.vcf

    bgzip -c diploidSV.vcf > ~{sample_id}.manta.vcf.gz
    tabix -p vcf ~{sample_id}.manta.vcf.gz
    
  >>>

  output {
    File vcf = "${sample_id}.manta.vcf.gz"
    File index = "${sample_id}.manta.vcf.gz.tbi"
  }

  runtime {
  cpus: "~{num_cpu_use}"
  requested_memory_mb_per_core: "~{mem_size_Mb}"
  runtime_minutes: "1800"
  }
#1400 minutes not enough for some samples
}

############
### Combine individual Manta VCFs
############
task CombineGvCFs {
  input {
    File GATK
    Array[File] input_gvcf
    Array[File] input_gvcf_index
    String output_vcf_filename
    File ref_fasta
    File ref_index
    File ref_dict
  }
  command <<<
    set -euo pipefail
    java -Xms8g \
      -jar ~{GATK} \
      CombineGVCFs \
      -R ~{ref_fasta} \
      -O "~{output_vcf_filename}.manta.gvcf" \
      ~{sep=' --variant ' input_gvcf }
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "8000"  
    runtime_minutes: "4700"  
  }

  output {
    File output_vcf = "~{output_vcf_filename}.manta.gvcf"
    File output_vcf_index = "~{output_vcf_filename}.manta.gvcf.tbi"
  }
}