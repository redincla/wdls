version 1.0

## Copyright CHUV, 2021
## Script to launch rare-variant tests

#################################################################
# WORKFLOW DEFINITION - RareVariantTest
#################################################################

workflow RareVariantTest {

String pipeline_version = "1.0"

input {
File input_vcf
File input_vcf_index
File GATK
File ref_fasta
File ref_index
File ref_dict
Int AN_threshold
Float AF_threshold
File gene_file
Array[String] single_test_names
Array[String] burden_test_names
File ped_file
}

String vcf_prefix=basename(input_vcf, ".vcf.gz")


### filter on observed AN: in 90% of called samples
call ANFilter {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      AN_threshold = AN_threshold,
      base_output_name = vcf_prefix
  }
  
### filter on AF: 0.01 and 0.05
call AFFilter {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      input_vcf = ANFilter.output_vcf,
      input_vcf_index = ANFilter.output_vcf_index,
      AF_threshold = AF_threshold,
      base_output_name = vcf_prefix
  }

### filter vcf for Himpact variants
    call HimpactFilter { 
      input:
        GATK = GATK,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_dict = ref_dict,
        input_vcf = AFFilter.output_vcf,
        input_vcf_index = AFFilter.output_vcf_index,
        base_output_name = vcf_prefix
      }

call convertForRVtest {
    input:
      input_vcf = HimpactFilter.output_vcf,
      input_vcf_index = HimpactFilter.output_vcf_index,
      prefix = vcf_prefix
}

scatter (idx in range(length(single_test_names))) {  
    call single_test {
    input:
      input_vcf = convertForRVtest.output_vcf,
      input_vcf_index = convertForRVtest.output_vcf_index,
      ped_file = ped_file,
      test_name = single_test_names[idx],
      base_output_name = vcf_prefix
  }
}

scatter (idx in range(length(burden_test_names))) {  
    call burden_test {
    input:
      input_vcf = convertForRVtest.output_vcf,
      input_vcf_index = convertForRVtest.output_vcf_index,
      gene_file = gene_file,
      ped_file = ped_file,
      test_name = burden_test_names[idx],
      base_output_name = vcf_prefix
  }
}

call convertForPLINK {
    input:
      input_vcf = HimpactFilter.output_vcf,
      input_vcf_index = HimpactFilter.output_vcf_index,
      ped_file = ped_file,
      prefix = vcf_prefix
}

call PLINKstd {
    input:
      input_bed = convertForPLINK.output_bed,
      input_bim = convertForPLINK.output_bim,
      input_fam = convertForPLINK.output_fam,
      ped_file = ped_file,
      prefix = vcf_prefix
}

call PLINKfisher {
    input:
      input_bed = convertForPLINK.output_bed,
      input_bim = convertForPLINK.output_bim,
      input_fam = convertForPLINK.output_fam,
      ped_file = ped_file,
      prefix = vcf_prefix
}


output {
    File filtered_vcf = HimpactFilter.output_vcf
    File filtered_vcf_index = HimpactFilter.output_vcf_index

    Array[File] single_test_assoc = single_test.output_assoc
    Array[File] single_test_logs = single_test.output_log

    Array[File] burden_test_assoc = burden_test.output_assoc
    Array[File] burden_test_logs = burden_test.output_log

    File PLINKstd_assoc = PLINKstd.output_assoc
    File PLINKstd_hh = PLINKstd.output_hh

    File PLINKfisher_fish = PLINKfisher.output_fish
}
}

  
#################################################################
# TASK DEFINITION 
#################################################################

############
### 1- filter on cohort AN : 90% of genotyped called
############
task ANFilter {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String base_output_name
    Int AN_threshold
  }

command <<<
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V ~{input_vcf} \
    -O ~{base_output_name}.highCR.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    -select "AN > ~{AN_threshold}"   
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "10000" 
  }

  output {
    File output_vcf = "~{base_output_name}.highCR.vcf.gz"
    File output_vcf_index = "~{base_output_name}.highCR.vcf.gz.tbi"
  }
}

############
###  2- filter on gnomAD AF
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
    gzip -cd ~{input_vcf} | sed 's@dbscSNV_ADA_SCORE,Number=.,Type=String@dbscSNV_ADA_SCORE,Number=1,Type=Float@g' | sed 's@dbscSNV_RF_SCORE,Number=.,Type=String@dbscSNV_RF_SCORE,Number=1,Type=Float@g' | sed 's@;dbscSNV_ADA_SCORE=\.;@;dbscSNV_ADA_SCORE=NaN;@g'  | sed 's@;dbscSNV_RF_SCORE=\.;@;dbscSNV_RF_SCORE=NaN;@g' | sed 's@;max_aaf_all=\.;@;max_aaf_all=NaN;@g' | bgzip -c > tmp.vcf.gz
    tabix tmp.vcf.gz
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V tmp.vcf.gz \
    -O ~{base_output_name}.highCR.lowAF.vcf.gz \
    -select "max_aaf_all < ~{AF_threshold} || max_aaf_all == 'NaN'"
>>> 

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.highCR.lowAF.vcf.gz"
    File output_vcf_index = "~{base_output_name}.highCR.lowAF.vcf.gz.tbi"
  }
}

############
### 3- Filter to retrieve High impact variants
############
task HimpactFilter {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String base_output_name
  }

command <<<
     java -Xmx8g -jar ~{GATK} \
     SelectVariants \
        -R ~{ref_fasta} \
        -V ~{input_vcf} \
        --keep-original-ac \
        -O "~{base_output_name}.highCR.lowAF.Himpact.vcf.gz" \
        -select "ExonicFunc.refGene == 'frameshift_deletion' || ExonicFunc.refGene == 'frameshift_insertion' || Func.refGene == 'splicing' || ExonicFunc.refGene == 'nonframeshift_deletion' || ExonicFunc.refGene == 'nonframeshift_insertion' || ExonicFunc.refGene == 'nonsynonymous_SNV' || ExonicFunc.refGene == 'startloss' || ExonicFunc.refGene == 'stopgain' || ExonicFunc.refGene == 'stoploss' || ExonicFunc.refGene == 'startgain' || dbscSNV_ADA_SCORE > 0.6 || dbscSNV_RF_SCORE > 0.6 || ClinVar_Sign  == 'Pathogenic/Likely_pathogenic' || ClinVar_Sign  == 'Likely_pathogenic' || ClinVar_Sign  == 'Pathogenic'"
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.highCR.lowAF.Himpact.vcf.gz"
    File output_vcf_index = "~{base_output_name}.highCR.lowAF.Himpact.vcf.gz.tbi"
  }
}

############
### 4- Convert vcf file for compatibility with rv-tests
############

task convertForRVtest {
  input {
    File input_vcf
    File input_vcf_index
    String prefix
  }

command <<<
  module add UHTS/Analysis/EPACTS/3.2.6
  gzip -cd ~{input_vcf} | sed  's@^chr@@g' | sed  's@ID=chr@ID=@g' | bgzip -c > ~{prefix}-rvtest.vcf.gz
  tabix ~{prefix}-rvtest.vcf.gz
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "2000" 
  runtime_minutes: "200"
  }

  output {
    File output_vcf = "~{prefix}-rvtest.vcf.gz"
    File output_vcf_index = "~{prefix}-rvtest.vcf.gz.tbi"
  }
}

############
### 5- launch single variant association test
############

task single_test {
  input {
    File input_vcf
    File input_vcf_index
    File ped_file
    String test_name
    String base_output_name
  }

command <<<
    /users/credin/refs/tools/rvtests/executable/rvtest --inVcf ~{input_vcf} \
    --pheno ~{ped_file} --out ~{base_output_name}_~{test_name} --single ~{test_name}
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  runtime_minutes: "200"
  }

  output {
    File output_assoc = "~{base_output_name}_~{test_name}.assoc"
    File output_log = "~{base_output_name}_~{test_name}.log"
  }
}

############
### 6- launch gene burden tests
############

task burden_test {
  input {
    File input_vcf
    File input_vcf_index
    File ped_file
    File gene_file
    String test_name
    String base_output_name
  }

command <<<
    /users/credin/refs/tools/rvtests/executable/rvtest --inVcf ~{input_vcf} \
    --pheno ~{ped_file} --out ~{base_output_name}_~{test_name} --geneFile ~{gene_file} \
    --burden ~{test_name} --vt price --freqUpper 0.05

>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "10000" 
  runtime_minutes: "400"
  }

  output {
    File output_assoc = "~{base_output_name}_~{test_name}.assoc"
    File output_log = "~{base_output_name}_~{test_name}.log"
  }
}

############
### 7- Convert vcf file for PLINK inputs
############

task convertForPLINK {
  input {
    File input_vcf
    File input_vcf_index
    File ped_file  ## FID should be 0, unless specific related individuals
    String prefix
  }

command <<<
  module add UHTS/Analysis/plink/1.90 
  cut -f 2 ~{ped_file} | sed  's@^@0\t@g' > sample.list
  plink --vcf ~{input_vcf} --const-fid --out ~{prefix}_PLINK --biallelic-only strict --keep-allele-order --indiv-sort file sample.list --make-bed
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "2000" 
  runtime_minutes: "200"
  }

  output {
    File output_bed = "~{prefix}_PLINK.bed"
    File output_bim = "~{prefix}_PLINK.bim"
    File output_fam = "~{prefix}_PLINK.fam"
  }
}


############
### 8- Launches PLINK stdandard association analyses
############

task PLINKstd {
  input {
    File input_bed
    File input_bim
    File input_fam
    File ped_file  ## FID should be 0, unless specific related individuals
    String prefix
  }

command <<<
  module add UHTS/Analysis/plink/1.90 
  plink --bed ~{input_bed} --bim ~{input_bim} --ci 0.95 --assoc --fam ~{ped_file} --out ~{prefix}_PLINK
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "2000" 
  runtime_minutes: "200"
  }

  output {
    File output_assoc = "~{prefix}_PLINK.assoc"
    File output_hh = "~{prefix}_PLINK.hh"
  }
}

############
### 9- Launches PLINK association analyses
############

task PLINKfisher {
  input {
    File input_bed
    File input_bim
    File input_fam
    File ped_file  ## FID should be 0, unless specific related individuals
    String prefix
  }

command <<<
  module add UHTS/Analysis/plink/1.90 
  plink --bed ~{input_bed} --bim ~{input_bim} --ci 0.95 --model --fisher --fam ~{ped_file} --out ~{prefix}_PLINK
>>>

  runtime {
  cpus: "1"
  requested_memory_mb_per_core: "2000" 
  runtime_minutes: "200"
  }

  output {
    File output_fish = "~{prefix}_PLINK.model"
  }
}