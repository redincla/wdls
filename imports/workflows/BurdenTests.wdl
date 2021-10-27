version 1.0

## Copyright CHUV, 2021
## Script to launch gene-burden tests
## Following workflow from https://github.com/DrGBL/WES.WGS

#################################################################
# WORKFLOW DEFINITION - BurdenTest
#################################################################

## Local import
import "/home/credin/scratch/WGS/wdls/imports/tasks/JointCalling-tasks-test.wdl" as Tasks

workflow BurdenTest {

String pipeline_version = "1.0"

input {
File input_vcf
File input_vcf_index
File ref_fasta
File ref_index
File ref_dict
File ped_file
File ancestry_IDs
String genome_version
File GATK
File tabix
File unpadded_intervals_file
Int scatter_count
Float AF_threshold
Int AC_threshold
File regions_list
}

String vcf_prefix=basename(input_vcf, ".vcf.gz")
Boolean is_hg38 = if ( genome_version == 'hg38' ) then true else false


  if (is_hg38) {
    call ConvertToGRCh38 {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        base_output_name = vcf_prefix
      }
  }

  call LeftAlign {
    input:
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      input_vcf = select_first([ConvertToGRCh38.output_vcf, input_vcf]),
      input_vcf_index = select_first([ConvertToGRCh38.output_vcf_index, input_vcf_index]),
      base_output_name = vcf_prefix
  }

  call HailQC {
    input:
      input_vcf = LeftAlign.output_vcf,
      input_vcf_index = LeftAlign.output_vcf_index,
      base_output_name = vcf_prefix
  }
  
  call AFFilter {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      input_vcf = HailQC.output_vcf,
      input_vcf_index = HailQC.output_vcf_index,
      AF_threshold = AF_threshold,
      base_output_name = vcf_prefix
  }

  call ACFilter {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      input_vcf = AFFilter.output_vcf,
      input_vcf_index = AFFilter.output_vcf_index,
      AC_threshold = AC_threshold,
      base_output_name = vcf_prefix
  }
  
  call SelectAncestry {
    input:
      ancestry_IDs = ancestry_IDs,
      input_vcf = ACFilter.output_vcf,
      input_vcf_index = ACFilter.output_vcf_index,
      base_output_name = vcf_prefix
  }

#  call PCACommon {
#    input:
#      input_vcf = SelectAncestry.output_vcf,
#      input_vcf_index = SelectAncestry.output_vcf_index,
#      base_output_name = vcf_prefix,
#      TMP_DIR = ''
#  }

#  call PCARare {
#    input:
#      input_vcf = SelectAncestry.output_vcf,
#      input_vcf_index = SelectAncestry.output_vcf_index,
#      base_output_name = vcf_prefix,
#      TMP_DIR = ''
#  }

  call MakeSitesOnlyVcf {
    input:
      input_vcf = SelectAncestry.output_vcf,
      input_vcf_index = SelectAncestry.output_vcf_index,
      base_output_name = vcf_prefix
  }

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
        input_vcf = MakeSitesOnlyVcf.output_vcf,
        input_vcf_index = MakeSitesOnlyVcf.output_vcf_index,
        interval = unpadded_intervals[idx],
        base_output_name = vcf_prefix,
        index = idx
    }

    call gvanno {
      input:
        input_vcf = SplitVCF.output_vcf,
        input_vcf_index = SplitVCF.output_vcf_index,
        base_output_name = vcf_prefix,
        genome_version = genome_version,
        index = idx,
        CURR_DIR = ''
    }

    call GetHighImpact as GetHighImpactALL {
      input:
        input_vcf = gvanno.output_vcf,
        input_vcf_index = gvanno.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }
  
    call GetModImpactNonMis as GetModImpactNonMisALL {
      input:
        input_vcf = gvanno.output_vcf,
        input_vcf_index = gvanno.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }

    call GetModImpactMissense as GetModImpactMissenseALL {
      input:
        input_vcf = gvanno.output_vcf,
        input_vcf_index = gvanno.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }

    call GetLowImpactSyn as GetLowImpactSynALL {
      input:
        input_vcf = gvanno.output_vcf,
        input_vcf_index = gvanno.output_vcf_index,
        base_output_name = vcf_prefix,
        index = idx
    }

    call GetHighImpact as GetHighImpactPASSONLY {
      input:
        input_vcf = gvanno.output_vcf_PASS,
        input_vcf_index = gvanno.output_vcf_index_PASS,
        base_output_name = vcf_prefix,
        index = idx
    }
  
    call GetModImpactNonMis as GetModImpactNonMisPASSONLY {
      input:
        input_vcf = gvanno.output_vcf_PASS,
        input_vcf_index = gvanno.output_vcf_index_PASS,
        base_output_name = vcf_prefix,
        index = idx
    }

    call GetModImpactMissense as GetModImpactMissensePASSONLY {
      input:
        input_vcf = gvanno.output_vcf_PASS,
        input_vcf_index = gvanno.output_vcf_index_PASS,
        base_output_name = vcf_prefix,
        index = idx
    }

    call GetLowImpactSyn as GetLowImpactSynPASSONLY {
      input:
        input_vcf = gvanno.output_vcf_PASS,
        input_vcf_index = gvanno.output_vcf_index_PASS,
        base_output_name = vcf_prefix,
        index = idx
    }
  }
  
  Array[File] HighImpactALL_input_vcfs = GetHighImpactALL.output_vcf
  Array[File] ModImpactNonMisALL_input_vcfs = GetModImpactNonMisALL.output_vcf
  Array[File] ModImpactMissenseALL_input_vcfs = GetModImpactMissenseALL.output_vcf
  Array[File] LowImpactSynALL_input_vcfs = GetLowImpactSynALL.output_vcf

  Array[File] HighImpactPASSONLY_input_vcfs = GetHighImpactPASSONLY.output_vcf
  Array[File] ModImpactNonMisPASSONLY_input_vcfs = GetModImpactNonMisPASSONLY.output_vcf
  Array[File] ModImpactMissensePASSONLY_input_vcfs = GetModImpactMissensePASSONLY.output_vcf
  Array[File] LowImpactSynPASSONLY_input_vcfs = GetLowImpactSynPASSONLY.output_vcf



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



  call Tasks.GatherVcfs as GatherHighImpactPASSONLY {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = HighImpactPASSONLY_input_vcfs,
      output_vcf_name = vcf_prefix + ".HighImpactPASSONLY.vcf.gz"
  }

  call Tasks.GatherVcfs as GatherModImpactNonMisPASSONLY {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = ModImpactNonMisPASSONLY_input_vcfs,
      output_vcf_name = vcf_prefix + ".ModImpactNonMisPASSONLY.vcf.gz"
  }

  call Tasks.GatherVcfs as GatherModImpactMissensePASSONLY {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = ModImpactMissensePASSONLY_input_vcfs,
      output_vcf_name = vcf_prefix + ".ModImpactMissensePASSONLY.vcf.gz"
  }

  call Tasks.GatherVcfs as GatherLowImpactSynPASSONLY {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = LowImpactSynPASSONLY_input_vcfs,
      output_vcf_name = vcf_prefix + ".LowImpactSynPASSONLY.vcf.gz"
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
        base_output_name = vcf_prefix + "LowImpactSynALL.qualifying.variants.list"
  }


  call CreateSNPFile as ListHighImpactPASSONLY {
    input:
        input_vcf = GatherHighImpactPASSONLY.output_vcf,
        input_vcf_index = GatherHighImpactPASSONLY.output_vcf_index,
        base_output_name = vcf_prefix + "HighImpactPASSONLY.qualifying.variants.list"
  }

  call CreateSNPFile as ListModImpactNonMisPASSONLY {
    input:
        input_vcf = GatherModImpactNonMisPASSONLY.output_vcf,
        input_vcf_index = GatherModImpactNonMisPASSONLY.output_vcf_index,
        base_output_name = vcf_prefix + "ModImpactNonMisPASSONLY.qualifying.variants.list"
  }

  call CreateSNPFile as ListModImpactMissensePASSONLY {
    input:
        input_vcf = GatherModImpactMissensePASSONLY.output_vcf,
        input_vcf_index = GatherModImpactMissensePASSONLY.output_vcf_index,
        base_output_name = vcf_prefix + "ModImpactMissensePASSONLY.qualifying.variants.list"
  }

  call CreateSNPFile as ListLowImpactSynPASSONLY {
    input:
        input_vcf = GatherLowImpactSynPASSONLY.output_vcf,
        input_vcf_index = GatherLowImpactSynPASSONLY.output_vcf_index,
        base_output_name = vcf_prefix + "LowImpactSynPASSONLY.qualifying.variants.list"
  }



  call CountCarriers as CountCarriersHighImpactPASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListHighImpactPASSONLY.variant_list,
        sample_list = sample_list,
        regions_list = regions_list,
        base_output_name = vcf_prefix + "HighImpactPASS.qualifying.variants.counts"
  }

  call CountCarriers as CountCarriersModImpactNonMisPASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListModImpactNonMisPASSONLY.variant_list,
        sample_list = sample_list,
        regions_list = regions_list,
        base_output_name = vcf_prefix + "ModImpactNonMisPASS.qualifying.variants.counts"
  }

    call CountCarriers as CountCarriersModImpactMissensePASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListModImpactMissensePASSONLY.variant_list,
        sample_list = sample_list,
        regions_list = regions_list,
        base_output_name = vcf_prefix + "ModImpactMissensePASS.qualifying.variants.counts"
  }

    call CountCarriers as CountCarriersLowImpactSynPASS {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variant_list = ListLowImpactSynPASSONLY.variant_list,
        sample_list = sample_list,
        regions_list = regions_list,
        base_output_name = vcf_prefix + "LowImpactSynPASS.calibrating.variants.counts"
  }

output {
  File hail_matrix = HailQC.hail_table
  File QC_vcf = HailQC.output_vcf
  File QC_vcf_index = HailQC.output_vcf_index

  File ancestry_vcf = SelectAncestry.output_vcf
  File ancestry_vcf_index = SelectAncestry.output_vcf_index

#  File common_eigenval = PCACommon.output_eigenval
#  File common_eigenvec = PCACommon.output_eigenvec
#  File common_log = PCACommon.output_log
#  File common_nosex = PCACommon.output_nosex

#  File rare_eigenval = PCARare.output_eigenval
#  File rare_eigenvec = PCARare.output_eigenvec
#  File rare_log = PCARare.output_log
#  File rare_nosex = PCARare.output_nosex

  File HighImpactPASSONLY_list = ListHighImpactPASSONLY.variant_list
  File ModImpactNonMisPASSONLY_list = ListModImpactNonMisPASSONLY.variant_list
  File ModImpactMissensePASSONLY_list = ListModImpactMissensePASSONLY.variant_list
  File LowImpactSynPASSONLY_list = ListLowImpactSynPASSONLY.variant_list

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
### 0- Convert from hg18 to GRCh38
############
task ConvertToGRCh38 {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
  }

command <<<
    module add UHTS/Analysis/EPACTS/3.2.6
    zgrep '^#' ~{input_vcf} | sed  's@contig=<ID=@contig=<ID=chr@g' > header
    gzip -cd ~{input_vcf} | grep -v '^#' | sed  's@^@chr@g' > ~{base_output_name}.GRCh38.vcf
    cat header ~{base_output_name}.GRCh38.vcf | bgzip -c > ~{base_output_name}.GRCh38.vcf.gz
    tabix -p vcf ~{base_output_name}.GRCh38.vcf.gz
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_vcf = "~{base_output_name}.GRCh38.vcf.gz"
    File output_vcf_index = "~{base_output_name}.GRCh38.vcf.gz.tbi"
  }
}


############
### 1- Left align variants
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
  cpus: "1"
	requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_vcf = "~{base_output_name}.normID.noChrM.vqsr.flt.vcf.gz"
    File output_vcf_index = "~{base_output_name}.normID.noChrM.vqsr.flt.vcf.gz.tbi"
  }
}

############
### 2- Format hail input table
############

task HailQC {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
  }

command <<<
    source /dcsrsoft/spack/bin/setup_dcsrsoft
    module load gcc python
    module load openblas/0.3.9
    module load mpich/3.3.2
    module load netlib-scalapack/2.1.0
    module load openjdk/1.8.0_222-b10

    TMP_DIR=`mktemp -d /tmp/tmp.XXXXXX`
    python3 <<CODE

    import hail as hl
    import hail.expr.aggregators as agg

    hl.init(spark_conf=None, tmp_dir= '${TMP_DIR}')
    hl.import_vcf('~{input_vcf}', min_partitions=4, reference_genome='GRCh38', force_bgz=True).write('~{base_output_name}.hail.full.normID.noChrM.mt', overwrite=True)
    metaData = hl.get_vcf_metadata('~{input_vcf}')

    mtAll = hl.read_matrix_table('~{base_output_name}.hail.full.normID.noChrM.mt')
    mtAll= mtAll.annotate_entries(AB = (mtAll.AD[1] / hl.sum(mtAll.AD) ))
    mtAll=hl.sample_qc(mtAll)
    mtAll = mtAll.filter_cols((mtAll.sample_qc.call_rate >= 0.97) & (mtAll.sample_qc.dp_stats.mean >= 20))
    mtAll = mtAll.filter_entries( (mtAll.GQ>=20) &
        (mtAll.DP >= 10) &
        ((mtAll.GT.is_hom_ref() & (mtAll.AB <= 0.1)) |
        (mtAll.GT.is_het() & (mtAll.AB >= 0.25) & (mtAll.AB <= 0.75)) |
        (mtAll.GT.is_hom_var() & (mtAll.AB >= 0.9))))
                        
    hl.export_vcf(mtAll, '~{base_output_name}.normID.GTflt.AB.noChrM.vqsr.flt.vcf.bgz', metadata=metaData, tabix=True)

    CODE
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File hail_table = "~{base_output_name}.hail.full.normID.noChrM.mt"
    File output_vcf = "~{base_output_name}.normID.GTflt.AB.noChrM.vqsr.flt.vcf.bgz"
    File output_vcf_index = "~{base_output_name}.normID.GTflt.AB.noChrM.vqsr.flt.vcf.bgz.tbi"
  }
}

############
### 3- Select Ancestry
############

task SelectAncestry {
  input {
    File ancestry_IDs
    File input_vcf
    File input_vcf_index
    String base_output_name
  }

command <<<
    module add UHTS/Analysis/samtools/1.10
    module add UHTS/Analysis/EPACTS/3.2.6

    bcftools view -S ~{ancestry_IDs} ~{input_vcf} -Oz > ~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.vqsr.flt.vcf.bgz
    tabix -p vcf ~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.vqsr.flt.vcf.bgz
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_vcf = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.vqsr.flt.vcf.bgz"
    File output_vcf_index = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.vqsr.flt.vcf.bgz.tbi"
  }
}

############
### 4a- PCA - Common variants
############

task PCACommon {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
    String TMP_DIR
  }

command <<<
    module add UHTS/Analysis/plink/1.90 
    TMP_DIR=`mktemp -d /tmp/tmp.XXXXXX`

    plink --vcf ~{input_vcf} \
    --biallelic-only strict \
    --chr 1-22 \
    --geno 0.05 \
    --snps-only 'just-acgt' \
    --hwe 1E-6 midp \
    --indep-pairwise 50 5 0.05 \
    --keep-allele-order \
    --mac 5 \
    --maf 0.01 \
    --out ${TMP_DIR}/commonAllelesPruned

    plink --vcf ~{input_vcf} \
    --extract ${TMP_DIR}/commonAllelesPruned.prune.in \
    --pca 10 \
    --out ~{base_output_name}.commonPCA

>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_eigenval = "~{base_output_name}.commonPCA.eigenval"
    File output_eigenvec = "~{base_output_name}.commonPCA.eigenvec"
    File output_log = "~{base_output_name}.commonPCA.log"
    File output_nosex = "~{base_output_name}.commonPCA.nosex"
  }
}

############
### 4b- PCA - Rare variants
############

task PCARare {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
    String TMP_DIR
  }

command <<<
    module add UHTS/Analysis/plink/1.90 
    TMP_DIR=`mktemp -d /tmp/tmp.XXXXXX`

    plink --vcf ~{input_vcf} \
    --biallelic-only strict \
    --chr 1-22 \
    --geno 0.05 \
    --snps-only 'just-acgt' \
    --hwe 1E-6 midp \
    --indep-pairwise 50 5 0.05 \
    --keep-allele-order \
    --max-maf 0.01 \
    --mac 2 \
    --out ${TMP_DIR}/rareAllelesPruned

    plink --vcf ~{input_vcf} \
    --extract ${TMP_DIR}/rareAllelesPruned.prune.in \
    --pca 20 \
    --out ~{base_output_name}.rarePCA
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_eigenval = "~{base_output_name}.rarePCA.eigenval"
    File output_eigenvec = "~{base_output_name}.rarePCA.eigenvec"
    File output_log = "~{base_output_name}.rarePCA.log"
    File output_nosex = "~{base_output_name}.rarePCA.nosex"
  }
}

############
###  5a- filter on gnomAD AF
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
    module add UHTS/Analysis/EPACTS/3.2.6
    gzip -cd ~{input_vcf} | sed 's@;max_aaf_all=\.;@;max_aaf_all=NaN;@g' | bgzip -c > tmp.vcf.gz
    tabix tmp.vcf.gz
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V tmp.vcf.gz \
    -O ~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAF.vqsr.flt.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    -select "max_aaf_all < ~{AF_threshold} || max_aaf_all == 'NaN'"
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
###  5b- filter on cohort AC
############
task ACFilter {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String base_output_name
    Float AC_threshold
  }

command <<<
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V ~{input_vcf} \
    -O ~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAFAC.vqsr.flt.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    -select "AC_Orig < ~{AC_threshold}"   
>>> 

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAFAC.vqsr.flt.vcf.gz"
    File output_vcf_index = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAFAC.vqsr.flt.vcf.gz.tbi"
  }
}

############
### 6 Make sites Only 
############

task MakeSitesOnlyVcf {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
  }
  command <<<
  module add UHTS/Analysis/EPACTS/3.2.6
  zcat ~{input_vcf} | grep "^#" | cut -f 1-8 > header
  zcat ~{input_vcf} | grep -v "^#" | cut -f 1-8 | awk 'BEGIN{FS=OFS="\t"} {gsub(".*", ".", $8)} 1' >  tmp.vcf
  cat header tmp.vcf | bgzip -c > ~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAFAC.vqsr.flt.sites-only.vcf.gz
  tabix -p vcf ~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAFAC.vqsr.flt.sites-only.vcf.gz
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "4000"   
    runtime_minutes: "60" 
  }

  output {
    File output_vcf = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAFAC.vqsr.flt.sites-only.vcf.gz"
    File output_vcf_index = "~{base_output_name}.normID.Ancestryflt.GTflt.AB.noChrM.lowAFAC.vqsr.flt.sites-only.vcf.gz.tbi"
  }
}


############
### 7- Split genome intervals for scattering
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
### 8- Split vcf
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
### 9- Annotation with VEP
############

task VEP {
  input {
    File input_vcf
    File input_vcf_index
    File interval ### insert this into VEP command!
    String base_output_name
  }

command <<<
    module add UHTS/Analysis/vep/96.0

    vep -i ~{input_vcf} \
    --plugin dbNSFP,VEP_canonical \
    -everything \
    --buffer_size 100000 \
    --force_overwrite \
    --offline \
    --fork 10 \
    --dir_cache $pathCache \
    --cache -o "~{base_output_name}".finalAnnot.txt

>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_file = "~{base_output_name}.finalAnnot.txt"
  }
}

############
### 9b- Annotation with VEP via gvanno
############
task gvanno {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
    String genome_version
    String CURR_DIR
    String index
  }

command <<<
  CURR_DIR=`pwd`
  source /dcsrsoft/spack/bin/setup_dcsrsoft
  module load gcc python
  source /home/credin/refs/tools/gvanno-1.3.2/venv/bin/activate
  cp ~{input_vcf} ${CURR_DIR}
  cp ~{input_vcf_index} ${CURR_DIR}
  cd ~/refs/tools/gvanno-1.3.2
  VCF=$(echo ~{input_vcf} | grep -o '[^/]*$')
  venv/bin/python ~/refs/tools/gvanno-1.3.2/gvanno.py ${CURR_DIR}/${VCF} ~/refs/tools/gvanno-1.3.2/ ${CURR_DIR} ~{genome_version} ~/refs/tools/gvanno-1.3.2/gvanno.toml ~{base_output_name}.~{index} --container singularity --no_vcf_validate --force_overwrite
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "12000" 
  }

  output {
    File output_vcf = "~{base_output_name}.~{index}_gvanno_grch38.vcf.gz"
    File output_vcf_index = "~{base_output_name}.~{index}_gvanno_grch38.vcf.gz.tbi"
    File output_vcf_PASS = "~{base_output_name}.~{index}_gvanno_pass_grch38.vcf.gz"
    File output_vcf_index_PASS = "~{base_output_name}.~{index}_gvanno_pass_grch38.vcf.gz.tbi"
  }
}


############
### 10a- Filter to retrieve High impact variants
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
    zgrep '^#' ~{input_vcf} | sed  's@contig=<ID=chr@contig=<ID=@g' > header
    zgrep -P '\|HIGH\|.[^,]*\|YES\|' ~{input_vcf} > ~{base_output_name}.HighImpact.~{index}.vcf
    cat header ~{base_output_name}.HighImpact.~{index}.vcf | bgzip -c > ~{base_output_name}.HighImpact.~{index}.vcf.gz
    tabix -p vcf ~{base_output_name}.HighImpact.~{index}.vcf.gz
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.HighImpact.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.HighImpact.~{index}.vcf.gz.tbi"
  }
}

############
### 10b- Filter to retrieve non missense Moderate impact variants
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
    zgrep '^#' ~{input_vcf} | sed  's@contig=<ID=chr@contig=<ID=@g' > header
    zgrep -E "inframe_insertion\|MODERATE\|.[^,]*\|YES\||inframe_deletion\|MODERATE\|.[^,]*\|YES\||protein_altering_variant\|MODERATE\|.[^,]*\|YES\|" ~{input_vcf} > ~{base_output_name}.ModImpactNonMis.~{index}.vcf
    cat header ~{base_output_name}.ModImpactNonMis.~{index}.vcf | bgzip -c > ~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz
    tabix -p vcf ~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.ModImpactNonMis.~{index}.vcf.gz.tbi"
  }
}

############
### 10c- Filter to retrieve Moderate impact - missense variants
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
    zgrep '^#' ~{input_vcf} | sed  's@contig=<ID=chr@contig=<ID=@g' > header
    zgrep -E "missense_variant\|MODERATE\|.[^,]*\|YES\|" ~{input_vcf} | grep -Ev "inframe_insertion|inframe_deletion|protein_altering" > ~{base_output_name}.ModImpactMis.~{index}.vcf
    cat header ~{base_output_name}.ModImpactMis.~{index}.vcf | bgzip -c > ~{base_output_name}.ModImpactMis.~{index}.vcf.gz
    tabix -p vcf ~{base_output_name}.ModImpactMis.~{index}.vcf.gz

>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.ModImpactMis.~{index}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.ModImpactMis.~{index}.vcf.gz.tbi"
  }
}


############
### 10d- Filter to retrieve synonymous variants - for calibration
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
    zgrep '^#' ~{input_vcf} | sed  's@contig=<ID=chr@contig=<ID=@g' > header
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
### 11- Gathers multiple VCF file
############
task GatherVcfs {
  input {
    File GATK
    File tabix
    Array[File] input_vcfs
    Array[File] input_vcfs_index
    String output_vcf_name
  }
    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    # check that output combined vcf is sorted, else tabix will fail
  command <<<
    set -euo pipefail
    module add UHTS/Analysis/samtools/1.10
    java -Xms6g \
      -jar ~{GATK} \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ~{sep=" --input " input_vcfs} \
      --output tmp_vcf.gz
    
    bcftools sort tmp_vcf.gz  -O z -o ~{output_vcf_name}

    ~{tabix} ~{output_vcf_name}
  >>>

  runtime {
    cpus: "1"
	requested_memory_mb_per_core: "7000" 
    runtime_minutes: "300"
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

############
### 12- Create SNP file for TRAPD
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

  python /users/credin/refs/tools/TRAPD/code/make_snp_file.py --vcffile ~{input_vcf} --vep --genecolname SYMBOL --outfile ~{base_output_name} 

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
### 13- Count carriers of qualifying variants
############

task CountCarriers {
  input {
    File input_vcf
    File input_vcf_index
    File variant_list
    File? sample_list
    String base_output_name
    File regions_list
  }

  String sample_list_arg = if sample_list then "--samplefile ~{sample_list}" else ""

command <<<
  source /dcsrsoft/spack/bin/setup_dcsrsoft
  module load gcc
  module load python/2.7.16

  python /users/credin/refs/tools/TRAPD/code/count_cases.py --vcffile ~{input_vcf} --snpfile ~{variant_list} --outfile ~{base_output_name} --pass ~{sample_list_arg} --bedfile ~{regions_list}

>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File count_file = "~{base_output_name}"
  }
}

############
### 14- Run Burden test
############

task RunBurdenTest {
  input {
    File cases_count_file
    File ctrls_count_file
    Int cases_size
    Int ctrls_size
    String base_output_name
  }

command <<<
  source /dcsrsoft/spack/bin/setup_dcsrsoft
  module load gcc
  module load r/4.0.2

  Rscript /users/credin/refs/tools/TRAPD/code/burden.R --casefile ~{cases_count_file} --casesize ~{cases_size} --controlfile ~{ctrls_count_file} --controlsize ~{ctrls_size} --outfile ~{base_output_name}
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
  File burden_file = "~{base_output_name}"
  }
}


############
### 8- Convert vcf file to Plink files
############

task VCFToPLINK {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
  }

command <<<
  module add UHTS/Analysis/plink/1.90 

  plink --vcf ~{input_vcf}  \
    --make-bed \
    --out "~{base_output_name}".Eur.normID.GTflt.AB.noChrM.vqsr.flt

  plink --bfile "~{base_output_name}".Eur.normID.GTflt.AB.noChrM.vqsr.flt \
    --hwe 1E-15 midp \
    --maf 0.01 \
    --geno 0.1 \
    --indep-pairwise 50 5 0.05 \
    --out "~{base_output_name}".regenie.LD.prune.maf1pct.geno10pct.Eur.normID.GTflt.AB.noChrM.vqsr.flt
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File binary_file = "~{base_output_name}.Eur.normID.GTflt.AB.noChrM.vqsr.flt"
    File Regenie_file = "~{base_output_name}.regenie.LD.prune.maf1pct.geno10pct.Eur.normID.GTflt.AB.noChrM.vqsr.flt"
  }
}