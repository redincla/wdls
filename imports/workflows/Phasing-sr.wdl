version 1.0

## Copyright CHUV, 2021
## Script to launch gene-burden tests
## Following workflow from https://github.com/DrGBL/WES.WGS

#################################################################
# WORKFLOW DEFINITION - Haplotype phasing from Illumina short reads
#################################################################

## Local import
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/JointCalling-tasks-test.wdl" as Tasks

workflow Phasing {

String pipeline_version = "1.0"

input {
    File sample_manifest # col 1: sample_id , col 2: path-to-single-sample-vcf , col3: path-to-single-sample-vcf-index, col4: path-to-cram, col 5: path-to-cram-index
    File chr_gmap #genetic recombination map /!\ careful to have the same format of chromosomes as in the input vcf: chr11 vs. 11
    File phased_ref  #bcf(.gz) csi.indexed file for phased reference genotypes  /!\ careful to have the same format of chromosomes as in the input vcf: chr11 vs. 11
    File phased_ref_index #bcf(.gz).csi indexed file
    String target_region #e.g. whole chromosome = chr20 ; specific range = chr20:2000000-3000000  /!\ careful to have the same format of chromosomes as in the input vcf: chr11 vs. 11
    File ref_fasta
    File ref_index
    File ref_dict
}
    Array[Array[String]] inputSamples = read_tsv(sample_manifest)

    scatter (line in inputSamples) {
        File sample_id = line[0]
        File input_vcf = line[1]
        String input_vcf_index = line[2]
        String input_cram = line[3]
        String input_cram_index = line[4]

    call whatshap_phasing {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        input_cram = input_cram,
        input_cram_index = input_cram_index,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        target_region = target_region,
        base_output_name = sample_id
    }

    call shapeit4 {
      input:
        input_vcf = whatshap_phasing.output_vcf,
        input_vcf_index = whatshap_phasing.output_vcf_index,
        base_output_name = sample_id,
        chr_gmap = chr_gmap,
        phased_ref = phased_ref,
        phased_ref_index = phased_ref_index,
        target_region = target_region
    }


    call get_gtf {
      input:
        shapeit4_phased_vcf = shapeit4.output_vcf,
        shapeit4_phased_vcf_index = shapeit4.output_vcf_index,
        whatshap_phased_vcf = whatshap_phasing.output_vcf,
        whatshap_phased_vcf_index = whatshap_phasing.output_vcf_index,
        base_output_name = sample_id
    }

    call phasing_reads as phasing_reads_shapeit4 {
      input:
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        input_cram = input_cram,
        input_cram_index = input_cram_index,
        phased_vcf = shapeit4.output_vcf,
        phased_vcf_index = shapeit4.output_vcf_index,
        base_output_name = sample_id + ".shapeit4"
    }

    call phasing_reads as phasing_reads_whatshap {
      input:
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        input_cram = input_cram,
        input_cram_index = input_cram_index,
        phased_vcf = whatshap_phasing.output_vcf,
        phased_vcf_index = whatshap_phasing.output_vcf_index,
        base_output_name = sample_id + ".whatshap"
    }

    call get_haplotypes as get_haplotypes_shapeit4 {
      input:
        target_region = target_region,
        phased_vcf = shapeit4.output_vcf,
        phased_vcf_index = shapeit4.output_vcf_index,
        base_output_name = sample_id
    }
}


output {
  Array[File] shapeit4_phased_vcfs = shapeit4.output_vcf
  Array[File] shapeit4_phased_vcfs_indices = shapeit4.output_vcf_index
  Array[File] shapeit4_phased_gtfs = get_gtf.output_gtf_shapeit4

  Array[File] whatshap_phased_vcfs = whatshap_phasing.output_vcf
  Array[File] whatshap_phased_vcfs_indices = whatshap_phasing.output_vcf_index
  Array[File] whatshap_phased_gtfs = get_gtf.output_gtf_whatshap

  Array[File] shapeit4_phased_bams = phasing_reads_shapeit4.phased_bam
  Array[File] shapeit4_phased_bams_indices = phasing_reads_shapeit4.phased_bam_index

  Array[File] whatshap_phased_bams = phasing_reads_whatshap.phased_bam
  Array[File] whatshap_phased_bams_indices = phasing_reads_whatshap.phased_bam_index

  Array[File] shapeit4_HapA = get_haplotypes_shapeit4.HapA
  Array[File] shapeit4_HapB = get_haplotypes_shapeit4.HapB
}
}


#################################################################
# TASK DEFINITION 
#################################################################

############
### 1- Haplotype phasing of target region using shapeit4
############

### /!\ If phasing a single sample, need to include a phased reference panel with --reference option
task shapeit4 {
  input {
    File input_vcf
    File input_vcf_index
    File chr_gmap  ## /!\ careful to have the same format of chromosomes as in the input vcf: chr11 vs. 11
    File phased_ref ## /!\ careful to have the same format of chromosomes as in the input vcf: chr11 vs. 11
    File phased_ref_index
    String target_region ## /!\ careful to have the same format of chromosomes as in the input vcf: chr11 vs. 11. whole chromosome = 20 ; specific range = 20:2000000-3000000
    String base_output_name
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load shapeit4/4.1.3
    module load htslib/1.12
    shapeit4 --input ~{input_vcf} --map ~{chr_gmap} --region ~{target_region} --thread 4 --sequencing --reference ~{phased_ref} --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --output "~{base_output_name}.shapeit4.phased.vcf.gz"
    tabix -p vcf ~{base_output_name}.shapeit4.phased.vcf.gz
>>>

  runtime {
    cpus: "2"
    requested_memory_mb_per_core: "30000" 
  }

  output {
    File output_vcf = "~{base_output_name}.shapeit4.phased.vcf.gz"
    File output_vcf_index = "~{base_output_name}.shapeit4.phased.vcf.gz.tbi"
  }
}

############
### 2- Haplotype phasing of target region using whatshap
############

task whatshap_phasing {
  input {
    File ref_fasta
    File ref_index
    File ref_dict
    String target_region ## /!\ careful to have the same format of chromosomes as in the input vcf: chr11 vs. 11
    File input_vcf
    File input_vcf_index
    File input_cram
    File input_cram_index
    String base_output_name
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load whatshap/1.0
    module load htslib/1.12
    whatshap phase -o "~{base_output_name}.whatshap.phased.vcf.gz" --indels --reference=~{ref_fasta} ~{input_vcf} ~{input_cram} --chromosome="~{target_region}"
    tabix -p vcf ~{base_output_name}.whatshap.phased.vcf.gz
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_vcf = "~{base_output_name}.whatshap.phased.vcf.gz"
    File output_vcf_index = "~{base_output_name}.whatshap.phased.vcf.gz.tbi"
  }
}

############
### 3- Getting gtf files from phased vcfs using whatshap
############

task get_gtf {
  input {
    File shapeit4_phased_vcf
    File shapeit4_phased_vcf_index
    File whatshap_phased_vcf
    File whatshap_phased_vcf_index
    String base_output_name
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load whatshap/1.0
    module load htslib/1.12

    whatshap stats --gtf="~{base_output_name}.shapeit4.phased.gtf" ~{shapeit4_phased_vcf}
    whatshap stats --gtf="~{base_output_name}.whatshap.phased.gtf" ~{whatshap_phased_vcf}
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_gtf_shapeit4 = "~{base_output_name}.shapeit4.phased.gtf"
    File output_gtf_whatshap = "~{base_output_name}.whatshap.phased.gtf"
  }
}


############
### 4- Phasing bam/cram files using whatshap
############

task phasing_reads {
  input {
    File ref_fasta
    File ref_index
    File ref_dict
    File input_cram
    File input_cram_index
    File phased_vcf
    File phased_vcf_index
    String base_output_name
  }

    Boolean is_bam = basename(input_cram, ".bam") + ".bam" == basename(input_cram)

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load whatshap/1.0
    module load htslib/1.12
    module load samtools/1.12

    if [[ "~{is_bam}" = true ]] 
      then
        whatshap haplotag -o "~{base_output_name}.chr11.haplotagged.bam" --reference=~{ref_fasta} ~{phased_vcf} ~{input_cram}
        samtools index "~{base_output_name}.chr11.haplotagged.bam"
      else
        samtools view ~{input_cram} -b -o "~{base_output_name}.chr11.bam" ## convert cram to bam on the fly  
        samtools index "~{base_output_name}.chr11.bam"
        whatshap haplotag -o "~{base_output_name}.chr11.haplotagged.bam" --reference=~{ref_fasta} ~{phased_vcf} "~{base_output_name}.chr11.bam"
        samtools index "~{base_output_name}.chr11.haplotagged.bam"
    fi  
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File phased_bam = "~{base_output_name}.chr11.haplotagged.bam"
    File phased_bam_index = "~{base_output_name}.chr11.haplotagged.bam.bai"
  }
}

############
### 5- Extract haplotypes 
############
task get_haplotypes {
  input {
    String target_region
    File phased_vcf
    File phased_vcf_index
    String base_output_name
  }

command <<<

    zgrep "^~{target_region}"  ~{phased_vcf}  | cut -f 1-5,10 | sed 's@|@\t@g' | cut -f 1-3,5-7 | cut -f 1-5 > ~{base_output_name}.shapeit4.HapA
    zgrep "^~{target_region}"  ~{phased_vcf}  | cut -f 1-5,10 | sed 's@|@\t@g' | cut -f 1-3,5-7 | cut -f 1-4,6 > ~{base_output_name}.shapeit4.HapB
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File HapA = "~{base_output_name}.shapeit4.HapA"
    File HapB = "~{base_output_name}.shapeit4.HapB"
  }
}
