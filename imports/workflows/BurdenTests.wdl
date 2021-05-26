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
  
  call SelectAncestry {
    input:
      ancestry_IDs = ancestry_IDs,
      input_vcf = HailQC.output_vcf,
      input_vcf_index = HailQC.output_vcf_index,
      base_output_name = vcf_prefix
  }

  call PCACommon {
    input:
      input_vcf = SelectAncestry.output_vcf,
      input_vcf_index = SelectAncestry.output_vcf_index,
      base_output_name = vcf_prefix,
      TMP_DIR = ''
  }

  call PCARare {
    input:
      input_vcf = SelectAncestry.output_vcf,
      input_vcf_index = SelectAncestry.output_vcf_index,
      base_output_name = vcf_prefix,
      TMP_DIR = ''
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
        input_vcf = SelectAncestry.output_vcf,
        input_vcf_index = SelectAncestry.output_vcf_index,
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
        TMP_DIR = ''
    }
  }

  Array[File] GatherVcfs_input_vcfs = gvanno.output_vcf

  call Tasks.GatherVcfs as GatherVcfs {
    input:
      GATK = GATK,
      tabix = tabix,
      input_vcfs = GatherVcfs_input_vcfs,
      output_vcf_name = vcf_prefix + "normID.Ancestryflt.GTflt.AB.noChrM.vqsr.flt.gvanno.vcf.bgz"
  }

output {
  File hail_matrix = HailQC.hail_table
  File QC_vcf = HailQC.output_vcf
  File QC_vcf_index = HailQC.output_vcf_index

  File ancestry_vcf = SelectAncestry.output_vcf
  File ancestry_vcf_index = SelectAncestry.output_vcf_index

  File common_PCAs = PCACommon.output_pcas
  File Rare_PCAs = PCARare.output_pcas

  File output_vcf_file = GatherVcfs.output_vcf
  File output_vcf_index_file = GatherVcfs.output_vcf_index
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
    --out ~{base_output_name}.commonPCA.txt

>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_pcas = "~{base_output_name}.commonPCA.txt"
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
    --out ~{base_output_name}.rarePCA.txt
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_pcas = "~{base_output_name}.rarePCA.txt"
  }
}

############
### 5- Split genome intervals for scattering
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
### 6- Split vcf
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
### 7- Annotation with VEP
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
    --plugin dbNSFP,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
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
### 7b- Annotation with VEP via gvanno
############
task gvanno {
  input {
    File input_vcf
    File input_vcf_index
    String base_output_name
    String genome_version
    String TMP_DIR
  }

command <<<
  TMP_DIR=`mktemp -d /tmp/tmp.XXXXXX`
  source /dcsrsoft/spack/bin/setup_dcsrsoft
  module load gcc python
  source /home/credin/refs/tools/gvanno-1.3.2/venv/bin/activate
  cd ~/refs/tools/gvanno-1.3.2
  venv/bin/python ~/refs/tools/gvanno-1.3.2/gvanno.py ~{input_vcf} ~/refs/tools/gvanno-1.3.2/ ${TMP_DIR} ~{genome_version} ~/refs/tools/gvanno-1.3.2/gvanno.toml ~{base_output_name} --container singularity --force_overwrite
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_vcf = "${TMP_DIR}/~{base_output_name}.finalAnnot.tsv.gz"
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

