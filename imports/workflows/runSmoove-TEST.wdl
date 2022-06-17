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

workflow smoove {
  # Run smoove SV detection algorithm on whole genomes, uses lumpy : https://github.com/brentp/smoove
  input {
    Array[File] input_bams #.bam file to search for SVs
    Array[File] input_bams_index #"[optional] Index for _cram_file. If not specified, index is assumed to be at cram_file_path + '.crai'
    Array[String] sample_ids #sample name. 
    File ref_fasta #.fasta file with reference used to align bam or cram file
    File ref_index #[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file
    File exclude_regions_bed #text file with lines specifying genomic intervals where SVs should not be called.Each line in (tab-separated) format
    String cohort_name
    File hg38_gff
  }

    scatter (idx in range(length(input_bams))) {  ## lumpy - used by smoove - only takes bam as inputs, to detect discordant readpairs

            call gSVCalling {  
                input:
                    ref_fasta = ref_fasta,
                    ref_index = ref_index,
                    exclude_regions_bed = exclude_regions_bed,
                    input_bam = input_bams[idx],
                    input_bam_index= input_bams_index[idx],
                    sample_id = sample_ids[idx]
            }

    File gSVCalls = gSVCalling.output_vcf
    File gSVCalls_index = gSVCalling.output_vcf_index
    }

    call getSVUnion {  
                input:
                    ref_fasta = ref_fasta,
                    ref_index = ref_index,
                    input_vcfs = gSVCalls,
                    input_vcfs_index = gSVCalls_index,
                    cohort_name = cohort_name
    }

    scatter (idx in range(length(input_bams))) {  
            call genotypeSV {
                input:
                    ref_fasta = ref_fasta,
                    ref_index = ref_index,
                    merged_vcf = getSVUnion.output_vcf,
                    input_bam = input_bams[idx],
                    input_bam_index = input_bams_index[idx],
                    sample_id = sample_ids[idx]
            }
    File genotypes = genotypeSV.output_vcf
    File genotypes_index = genotypeSV.output_vcf_index
    }

    call MergeVCF {
            input:
                genotyped_vcfs = genotypes,
                genotyped_index = genotypes_index,
                cohort_name = cohort_name
    }


    call AnnotateVCF  {
            input:
                hg38_gff = hg38_gff,
                genotyped_vcf = MergeVCF.output_vcf,
                genotyped_vcf_index = MergeVCF.output_vcf_index,
                cohort_name = cohort_name
    }

  output {
    File final_vcf = MergeVCF.output_vcf
    File final_vcf_index =  MergeVCF.output_vcf_index
    File final_vcf_annotated = AnnotateVCF.output_vcf
    File final_vcf_annotated_index =  AnnotateVCF.output_vcf_index
  }
}

############
### 2- Launch smoove germline SV caller on individual samples
############

task gSVCalling {
  input {
    File ref_fasta
    File ref_index
    File exclude_regions_bed
    File input_bam
    File input_bam_index
    String sample_id   
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load singularity
    export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/data"

    singularity run /dcsrsoft/singularity/containers/smoove-0.2.7.sif smoove \
    call --outdir ./ --exclude ~{exclude_regions_bed} \
    --name ~{sample_id} --fasta ~{ref_fasta} \
    -p 1 --genotype ~{input_bam}
  >>>

  output {
    File output_vcf = "~{sample_id}-smoove.genotyped.vcf.gz"
    File output_vcf_index = "~{sample_id}-smoove.genotyped.vcf.gz.csi"
  }

  runtime {
    cpus: "1"
    requested_memory_mb_per_core: "10000"  #10000 for large WGS initially ~{mem_size_Mb}, trying to increase for stuck samples
    runtime_minutes: "400"
  }
}

############
### 3- Get the union of SV sites across all samples
############

task getSVUnion {
  input {
    File ref_fasta
    File ref_index
    Array[File] input_vcfs
    Array[File] input_vcfs_index
    String cohort_name   
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load singularity
    export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/data"

    vcfInputsList=~{write_lines(input_vcfs)}
    vcfIndexInputsList=~{write_lines(input_vcfs_index)}

    cp $vcfInputsList vcfInputs.list
    cp $vcfIndexInputsList vcfIndexInputs.list

    while read line; do
     SAMPLEID=$(basename ${line} .genotyped.vcf.gz)
     ln -s $line "${SAMPLEID}.genotyped.vcf.gz"
    done < vcfInputs.list
    singularity run /dcsrsoft/singularity/containers/smoove-0.2.7.sif smoove \
    merge --name ~{cohort_name} -f ~{ref_fasta} --outdir ./ ./*.genotyped.vcf.gz   
>>>

  output {
    File output_vcf = "~{cohort_name}.sites.vcf.gz"
  }

  runtime {
    cpus: "1"
    requested_memory_mb_per_core: "10000"
    runtime_minutes: "400"
  }
}

# if gathered array file does not work as input file:create symbolic link for all files in workind directory
# while read line; do
#    root=$(basename ~{vcf_file} .vcf.gz)
#    ln -s ~{vcf_file} "~{root}.vcf.gz
# done < gvcf


############
### 4- Genotype each individual sample using the cohort results as backbone
############

task genotypeSV {
  input {
    File ref_fasta
    File ref_index
    File merged_vcf
    File input_bam
    File input_bam_index
    String sample_id   
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load singularity
    export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/data"

    singularity run /dcsrsoft/singularity/containers/smoove-0.2.7.sif smoove \
    genotype -d -x -p 1 --name ~{sample_id} \
    --outdir ./ --fasta ~{ref_fasta} \
    --vcf ~{merged_vcf} ~{input_bam}
>>>

  output {
    File output_vcf = "~{sample_id}-smoove.genotyped.vcf.gz"
    File output_vcf_index = "~{sample_id}-smoove.genotyped.vcf.gz.csi"
  }

 runtime {
    cpus: "1"
    requested_memory_mb_per_core: "10000"
    runtime_minutes: "400"
  }
}

############
### 5- Generate combined gVCF for the cohort
############
task MergeVCF {
  input {
    Array[File] genotyped_vcfs
    Array[File] genotyped_index
    String cohort_name   
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load singularity
    module load htslib/1.12
    export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/data"

    vcfInputsList=~{write_lines(genotyped_vcfs)}
    vcfIndexInputsList=~{write_lines(genotyped_index)}

    cp $vcfInputsList vcfInputs.list
    cp $vcfIndexInputsList vcfIndexInputs.list

    while read line; do
     SAMPLEID=$(basename ${line} .genotyped.vcf.gz)
     ln -s $line "${SAMPLEID}.genotyped.vcf.gz"
     ln -s ${line}.csi "${SAMPLEID}.genotyped.vcf.gz.csi"
    done < vcfInputs.list

    singularity run /dcsrsoft/singularity/containers/smoove-0.2.7.sif smoove \
    paste --name ~{cohort_name} ./*.vcf.gz
    tabix "~{cohort_name}.smoove.square.vcf.gz"
>>>

  output {
    File output_vcf = "~{cohort_name}.smoove.square.vcf.gz"
    File output_vcf_index = "~{cohort_name}.smoove.square.vcf.gz.tbi"
  }

 runtime {
    cpus: "1"
    requested_memory_mb_per_core: "10000"
    runtime_minutes: "400"
  }
}

############
### 6- Annotate SVs (preliminary annotations only)
############

task AnnotateVCF {
  input {
    File genotyped_vcf
    File genotyped_vcf_index
    File hg38_gff
    String cohort_name   
  }

command <<<    
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc
    module load singularity
    module load htslib/1.12
    export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/data"
    
    cp ~{genotyped_vcf} "./~{cohort_name}.smoove.square.vcf.gz"

    singularity run /dcsrsoft/singularity/containers/smoove-0.2.7.sif smoove \
    annotate --gff ~{hg38_gff} ~{cohort_name}.smoove.square.vcf.gz | bgzip -c > "~{cohort_name}.smoove.square.anno.vcf.gz"
    tabix "~{cohort_name}.smoove.square.anno.vcf.gz"
>>>

  output {
    File output_vcf = "~{cohort_name}.smoove.square.anno.vcf.gz"
    File output_vcf_index = "~{cohort_name}.smoove.square.anno.vcf.tbi"
  }

 runtime {
    cpus: "1"
    requested_memory_mb_per_core: "10000"
    runtime_minutes: "400"
  }
}