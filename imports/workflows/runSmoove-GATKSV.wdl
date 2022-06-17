##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-sv/Delly/1/wdl
## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements SV calling with DELLY2 by Tobias Rausch (https://github.com/dellytools/delly)
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

workflow Smoove {
  # Run smoove SV detection algorithm on whole genomes, uses lumpy : https://github.com/brentp/smoove
  input {
    File bam_or_cram_file #.bam file to search for SVs
    File bam_or_cram_index #"[optional] Index for _cram_file. If not specified, index is assumed to be at cram_file_path + '.crai'
    String sample_id #sample name. 
    File reference_fasta #.fasta file with reference used to align bam or cram file
    File reference_index #[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file
    File region_bed #text file with lines specifying genomic intervals where SVs should not be called.Each line in (tab-separated) format
  }

            call gSVCalling {  
                input:
                    reference_fasta = reference_fasta,
                    reference_index = reference_index,
                    region_bed = region_bed,
                    bam_or_cram_file = bam_or_cram_file,
                    bam_or_cram_index = bam_or_cram_index,
                    sample_id = sample_id
            }

  output {
    File vcf = gSVCalling.output_vcf
    File index = gSVCalling.output_vcf_index
  }
}

############
### Launch smoove germline SV caller on individual samples
############

task gSVCalling {
  input {
    File reference_fasta
    File reference_index
    File region_bed
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id   
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load singularity
    export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/data"

    singularity run /dcsrsoft/singularity/containers/smoove-0.2.7.sif smoove \
    call --outdir ./ --exclude ~{region_bed} \
    --name ~{sample_id} --fasta ~{reference_fasta} \
    -p 1 --genotype ~{bam_or_cram_file}
  >>>

  output {
    File output_vcf = "~{sample_id}-smoove.genotyped.vcf.gz"
    File output_vcf_index = "~{sample_id}-smoove.genotyped.vcf.gz.csi"
  }

  runtime {
    cpus: "1"
    requested_memory_mb_per_core: "5000"  #10000 for large WGS initially ~{mem_size_Mb}, trying to increase for stuck samples
    runtime_minutes: "120"
  }
}