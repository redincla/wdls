version 1.0

## Copyright CHUV, 2021
## Script to prepare for gene-burden tests using population cohort
## Following workflow from https://github.com/DrGBL/WES.WGS and TRAPD https://github.com/mhguo1/TRAPD

#################################################################
# WORKFLOW DEFINITION - BurdenTest
#################################################################

## Local import
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/JointCalling-tasks-test.wdl" as Tasks

workflow BurdenTestgnomAD {

String pipeline_version = "1.0"

input {
File input_vcf
File input_vcf_index
File ref_fasta
File ref_index
File ref_dict
String genome_version
File GATK
File tabix
File unpadded_intervals_file
Int scatter_count
Float AF_threshold
File regions_list
}

String vcf_prefix=basename(input_vcf, ".vcf.gz")

call SplitIntervalList {
    input:
      GATK = GATK,
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict
  }

Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    call SplitVCF {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        interval = unpadded_intervals[idx],
        base_output_name = vcf_prefix,
        index = idx
    }
  
    call LeftAlign {
    input:
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      input_vcf = SplitVCF.output_vcf,
      input_vcf_index = SplitVCF.output_vcf_index,
      base_output_name = vcf_prefix
  }
  
  call AFFilter {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      input_vcf = LeftAlign.output_vcf,
      input_vcf_index = LeftAlign.output_vcf_index,
      AF_threshold = AF_threshold,
      base_output_name = vcf_prefix
  }

    call GetHighImpact as GetHighImpactALL {
      input:
        input_vcf = AFFilter.output_vcf,
        input_vcf_index = AFFilter.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }
  
    call GetModImpactNonMis as GetModImpactNonMisALL {
      input:
        input_vcf = AFFilter.output_vcf,
        input_vcf_index = AFFilter.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }

    call GetModImpactMissense as GetModImpactMissenseALL {
      input:
        input_vcf = AFFilter.output_vcf,
        input_vcf_index = AFFilter.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }

    call GetLowImpactSyn as GetLowImpactSynALL {
      input:
        input_vcf = AFFilter.output_vcf,
        input_vcf_index = AFFilter.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }
  }
  
  Array[File] HighImpactALL_input_vcfs = GetHighImpactALL.output_vcf
  Array[File] ModImpactNonMisALL_input_vcfs = GetModImpactNonMisALL.output_vcf
  Array[File] ModImpactMissenseALL_input_vcfs = GetModImpactMissenseALL.output_vcf
  Array[File] LowImpactSynALL_input_vcfs = GetLowImpactSynALL.output_vcf

  call Tasks.GatherVcfs as GatherHighImpactALL {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = HighImpactALL_input_vcfs,
      output_vcf_name = vcf_prefix + ".HighImpactALL.vcf.gz"
  }

  call Tasks.GatherVcfs as GatherModImpactNonMisALL {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = ModImpactNonMisALL_input_vcfs,
      output_vcf_name = vcf_prefix + ".ModImpactNonMisALL.vcf.gz"
  }

  call Tasks.GatherVcfs as GatherModImpactMissenseALL {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = ModImpactMissenseALL_input_vcfs,
      output_vcf_name = vcf_prefix + ".ModImpactMissenseALL.vcf.gz"
  }

  call Tasks.GatherVcfs as GatherLowImpactSynALL {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = LowImpactSynALL_input_vcfs,
      output_vcf_name = vcf_prefix + ".LowImpactSynALL.vcf.gz"
  }

  call CreateSNPFile as ListHighImpactALL {
    input:
        input_vcf = GatherHighImpactALL.output_vcf,
        input_vcf_index = GatherHighImpactALL.output_vcf_index,
        base_output_name = vcf_prefix + "HighImpactALL.qualifying.variants.list"
  }

    call CreateSNPFile as ListModImpactNonMisALL {
    input:
        input_vcf = GatherModImpactNonMisALL.output_vcf,
        input_vcf_index = GatherModImpactNonMisALL.output_vcf_index,
        base_output_name = vcf_prefix + "ModImpactNonMisALL.qualifying.variants.list"
  }

    call CreateSNPFile as ListModImpactMissenseALL {
    input:
        input_vcf = GatherModImpactMissenseALL.output_vcf,
        input_vcf_index = GatherModImpactMissenseALL.output_vcf_index,
        base_output_name = vcf_prefix + "ModImpactMissenseALL.qualifying.variants.list"
  }

    call CreateSNPFile as ListLowImpactSynALL {
    input:
        input_vcf = GatherLowImpactSynALL.output_vcf,
        input_vcf_index = GatherLowImpactSynALL.output_vcf_index,
        base_output_name = vcf_prefix + "LowImpactSynALL.calibrating.variants.list"
  }

  call CountCarriers as CountCarriersHighImpactPASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListHighImpactALL.variant_list,
        base_output_name = vcf_prefix + "HighImpactPASS.qualifying.variants.counts"
  }

  call CountCarriers as CountCarriersModImpactNonMisPASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListModImpactNonMisALL.variant_list,
        base_output_name = vcf_prefix + "ModImpactNonMisPASS.qualifying.variants.counts"
  }

    call CountCarriers as CountCarriersModImpactMissensePASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListModImpactMissenseALL.variant_list,
        base_output_name = vcf_prefix + "ModImpactMissensePASS.qualifying.variants.counts"
  }

    call CountCarriers as CountCarriersLowImpactSynPASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListLowImpactSynALL.variant_list,
        base_output_name = vcf_prefix + "LowImpactSyn.calibrating.variants.counts"
  }


output {
  File HighImpactALL_list = ListHighImpactALL.variant_list
  File ModImpactNonMisALL_list = ListModImpactNonMisALL.variant_list
  File ModImpactMissenseALL_list = ListModImpactMissenseALL.variant_list
  File LowImpactSynALL_list = ListLowImpactSynALL.variant_list

  File HighImpactPASS_count = CountCarriersHighImpactPASS.count_file
  File ModImpactNonMisPASS_count = CountCarriersModImpactNonMisPASS.count_file
  File ModImpactMissensePASS_count = CountCarriersModImpactMissensePASS.count_file
  File LowImpactSynPASS_count = CountCarriersLowImpactSynPASS.count_file
} 
}


#################################################################
# TASK DEFINITION 
#################################################################

############
### 1- Split genome intervals for scattering
############
task SplitIntervalList {
  input {
    File GATK
    File interval_list
    Int scatter_count
    File ref_fasta
    File ref_index
    File ref_dict
    String scatter_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
  }
  command {
    java -Xms3g \
      -jar ~{GATK} \
      SplitIntervals \
      -L ~{interval_list} -O scatterDir -scatter ~{scatter_count} -R ~{ref_fasta} \
      -imr OVERLAPPING_ONLY \
      -mode ~{scatter_mode}
   }

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "3500"
    runtime_minutes: "20"
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

############
### 2- Split vcf
############
task SplitVCF {
  input {
    File input_vcf
    File input_vcf_index
    File interval
    String base_output_name
    String index
  }

  command <<<
  module add UHTS/Analysis/samtools/1.10
  module add UHTS/Analysis/EPACTS/3.2.6
  grep '^chr' ~{interval} > "~{interval}_clean"
  bcftools view -R "~{interval}_clean" ~{input_vcf} -O z -o "~{base_output_name}.GRCh38.~{index}.vcf.gz"
  tabix -p vcf "~{base_output_name}.GRCh38.~{index}.vcf.gz"
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "6000"
    runtime_minutes: "200"
  }

  output {
    File output_vcf = "~{base_output_name}.GRCh38.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.GRCh38.~{index}.vcf.gz.tbi"
  }
}

############
### 3- Left align variants
############
task LeftAlign {
  input {
    File ref_fasta
    File ref_index
    File input_vcf
    File input_vcf_index
    String base_output_name
  }

command <<<
    module add UHTS/Analysis/samtools/1.10
    module add UHTS/Analysis/EPACTS/3.2.6
    bcftools filter -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX \
        ~{input_vcf} -Ou | \
        bcftools norm -m -any --check-ref w -f ~{ref_fasta} -Ou | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > "~{base_output_name}.normID.noChrM.vqsr.flt.vcf.gz"
    tabix -p vcf ~{base_output_name}.normID.noChrM.vqsr.flt.vcf.gz
>>>

  runtime {
  runtime_minutes: "200"
  cpus: "1"
	requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_vcf = "~{base_output_name}.normID.noChrM.vqsr.flt.vcf.gz"
    File output_vcf_index = "~{base_output_name}.normID.noChrM.vqsr.flt.vcf.gz.tbi"
  }
}

############
###  4- filter on gnomAD AF
############
task AFFilter {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String base_output_name
    Float AF_threshold
  }

command <<<
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V ~{input_vcf} \
    -O ~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAF.vqsr.flt.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    -select "AF < ~{AF_threshold}"
>>> 

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAF.vqsr.flt.vcf.gz"
    File output_vcf_index = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAF.vqsr.flt.vcf.gz.tbi"
  }
}


############
### 5a- Filter to retrieve High impact variants
############
task GetHighImpact {
  input {
    File input_vcf
    File input_vcf_index
    String index
    String base_output_name
  }

# NULL_VARIANT=True, LoF=HC,  , 
command <<<
    module add UHTS/Analysis/EPACTS/3.2.6
    zgrep '^#' ~{input_vcf} > header
    zgrep -P '\|HIGH\|.[^,]*\|YES\|' ~{input_vcf} > ~{base_output_name}.HighImpact.~{index}.vcf
    cat header ~{base_output_name}.HighImpact.~{index}.vcf | bgzip -c > ~{base_output_name}.HighImpact.~{index}.vcf.gz
    tabix -p vcf ~{base_output_name}.HighImpact.~{index}.vcf.gz
>>>

  runtime {
  runtime_minutes: "200"  
  cpus: "1"
	requested_memory_mb_per_core: "12000" 
  }

  output {
    File output_vcf = "~{base_output_name}.HighImpact.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.HighImpact.~{index}.vcf.gz.tbi"
  }
}

############
### 5b- Filter to retrieve non missense Moderate impact variants
############
task GetModImpactNonMis {
  input {
    File input_vcf
    File input_vcf_index
    String index
    String base_output_name
  }

command <<<
    module add UHTS/Analysis/EPACTS/3.2.6
    zgrep '^#' ~{input_vcf} > header
    zgrep -E "inframe_insertion\|MODERATE\|.[^,]*\|YES\||inframe_deletion\|MODERATE\|.[^,]*\|YES\||protein_altering_variant\|MODERATE\|.[^,]*\|YES\|" ~{input_vcf} > ~{base_output_name}.ModImpactNonMis.~{index}.vcf
    cat header ~{base_output_name}.ModImpactNonMis.~{index}.vcf | bgzip -c > ~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz
    tabix -p vcf ~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "12000" 
  }

  output {
    File output_vcf = "~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz.tbi"
  }
}

############
### 5c- Filter to retrieve Moderate impact - missense variants
############
task GetModImpactMissense {
  input {
    File input_vcf
    File input_vcf_index
    String index
    String base_output_name
  }

command <<<
    module add UHTS/Analysis/EPACTS/3.2.6
    zgrep '^#' ~{input_vcf} > header
    zgrep -E "missense_variant\|MODERATE\|.[^,]*\|YES\|" ~{input_vcf} | grep -Ev "inframe_insertion|inframe_deletion|protein_altering" > ~{base_output_name}.ModImpactMis.~{index}.vcf
    cat header ~{base_output_name}.ModImpactMis.~{index}.vcf | bgzip -c > ~{base_output_name}.ModImpactMis.~{index}.vcf.gz
    tabix -p vcf ~{base_output_name}.ModImpactMis.~{index}.vcf.gz

>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "12000" 
  }

  output {
    File output_vcf = "~{base_output_name}.ModImpactMis.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.ModImpactMis.~{index}.vcf.gz.tbi"
  }
}

############
### 5d- Filter to retrieve synonymous variants - for calibration
############
task GetLowImpactSyn {
  input {
    File input_vcf
    File input_vcf_index
    String index
    String base_output_name
  }

command <<<
    module add UHTS/Analysis/EPACTS/3.2.6
    zgrep '^#' ~{input_vcf}  > header
    zgrep -E "synonymous_variant\|LOW\|.[^,]*\|YES\|" ~{input_vcf} | grep -Ev "inframe_insertion|inframe_deletion|protein_altering|missense_variant" > ~{base_output_name}.LowImpactSyn.~{index}.vcf
    cat header ~{base_output_name}.LowImpactSyn.~{index}.vcf | bgzip -c > ~{base_output_name}.LowImpactSyn.~{index}.vcf.gz
    tabix -p vcf ~{base_output_name}.LowImpactSyn.~{index}.vcf.gz

>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.LowImpactSyn.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.LowImpactSyn.~{index}.vcf.gz.tbi"
  }
}

############
### 6- Create SNP file for TRAPD
############

task CreateSNPFile {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
  }

command <<<
  source /dcsrsoft/spack/bin/setup_dcsrsoft
  module load gcc
  module load python/2.7.16

  python /scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/redin/tools/TRAPD/code/make_snp_file.py --vcffile ~{input_vcf} --vep --genecolname SYMBOL --outfile ~{base_output_name}

>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File variant_list = "~{base_output_name}"
  }
}

############
### 7- Count carriers of qualifying variants
############

task CountCarriers {
  input {
    File input_vcf
    File input_vcf_index
    File variant_list
    String base_output_name
    File regions_list
  }

command <<<
  source /dcsrsoft/spack/bin/setup_dcsrsoft
  module load gcc
  module load python/2.7.16

  python /scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/redin/tools/TRAPD/code/count_controls.py --vcffile ~{input_vcf} --snpfile ~{variant_list} --outfile ~{base_output_name} --pass --database gnomad --bedfile ~{regions_list}
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "12000" 
  }

  output {
    File count_file = "~{base_output_name}"
  }
}
