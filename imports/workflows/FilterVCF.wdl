version 1.0

## Copyright CHUV, 2020
## Script to split multisample vcf files 
## and perform automatic pre-filter

#################################################################
# WORKFLOW DEFINITION - VCF Filtering
#################################################################

workflow VCFSplitAndFilter {

String pipeline_version = "1.0"

input {
Array[String] sample_list
File input_vcf
File input_vcf_index
File GATK
File ref_fasta
File ref_index
File ref_dict
}

String base_name = basename(input_vcf, ".vcf.gz")
#Array[String] sample_list = read_tsv(sample_IDs)

call ACFilter {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      base_output_name = base_name
  }

call AFFilter {
    input:
      GATK = GATK,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_dict = ref_dict,
      input_vcf = ACFilter.output_vcf,
      input_vcf_index = ACFilter.output_vcf_index,
      base_output_name = base_name
  }

scatter (idx in range(length(sample_list))) {  
    call SplitAndHQFilter { 
        input:
            GATK = GATK,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            input_vcf = AFFilter.output_vcf,
            input_vcf_index = AFFilter.output_vcf_index,
            sample_ID = sample_list[idx],
            base_output_name = base_name
        }
    
    call HimpactFilter { 
        input:
            GATK = GATK,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            input_vcf = SplitAndHQFilter.output_vcf,
            input_vcf_index = SplitAndHQFilter.output_vcf_index,
            sample_ID = sample_list[idx],
            base_output_name = base_name
        }
    
    call VCFToTSV { 
        input:
            GATK = GATK,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            input_vcf = HimpactFilter.output_vcf,
            input_vcf_index = HimpactFilter.output_vcf_index,
            sample_ID = sample_list[idx],
            base_output_name = base_name
    }
}

output {
    Array[File] scattered_vcfs = HimpactFilter.output_vcf
    Array[File] scattered_vcfs_index = HimpactFilter.output_vcf_index
    Array[File] scattered_tables = VCFToTSV.output_tsv
  }

}

#################################################################
# TASK DEFINITION 
#################################################################


############
### filter on cohort AC
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
  }

command <<<
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V ~{input_vcf} \
    -O ~{base_output_name}.lowAC.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    -select "AC < 6"   
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.lowAC.vcf.gz"
    File output_vcf_index = "~{base_output_name}.lowAC.vcf.gz.tbi"
  }
}

############
### filter on databases AF
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
  }

command <<<
    module add UHTS/Analysis/EPACTS/3.2.6
    gzip -cd ~{input_vcf} | sed 's@AF_popmax,Number=.,Type=String@AF_popmax,Number=1,Type=Float@g' | sed 's@dbscSNV_ADA_SCORE,Number=.,Type=String@dbscSNV_ADA_SCORE,Number=1,Type=Float@g' | sed 's@dbscSNV_RF_SCORE,Number=.,Type=String@dbscSNV_RF_SCORE,Number=1,Type=Float@g' | sed 's@dpsi_max_tissue,Number=.,Type=String@dpsi_max_tissue,Number=1,Type=Float@g' | sed 's@dpsi_zscore,Number=.,Type=String@dpsi_zscore,Number=1,Type=Float@g' | sed 's@;AF_popmax=\.;@;AF_popmax=NaN;@g' | sed 's@;dbscSNV_ADA_SCORE=\.;@;dbscSNV_ADA_SCORE=NaN;@g'  | sed 's@;dbscSNV_RF_SCORE=\.;@;dbscSNV_RF_SCORE=NaN;@g' | sed 's@;dpsi_max_tissue=\.;@;dpsi_max_tissue=NaN;@g' | sed 's@;dpsi_zscore=\.;@;dpsi_zscore=NaN;@g' | bgzip -c > tmp.vcf.gz
    tabix tmp.vcf.gz
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V tmp.vcf.gz \
    -O ~{base_output_name}.lowAC.lowAF.vcf.gz \
    -select "AF_popmax < 0.01 || AF_popmax == 'NaN'"
>>>  # need to push everything under a single expression with || else only takes into account the last select expression??

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.lowAC.lowAF.vcf.gz"
    File output_vcf_index = "~{base_output_name}.lowAC.lowAF.vcf.gz.tbi"
  }
}

############
### Split by sample && filter on GQ/DP
############
task SplitAndHQFilter {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String sample_ID
    String base_output_name
  }

command <<<
     java -Xmx8g -jar ~{GATK} \
     SelectVariants \
        -R ~{ref_fasta} \
        -V ~{input_vcf} \
        -O "~{base_output_name}.lowAC.lowAF.HQ.~{sample_ID}.vcf.gz" \
        --exclude-non-variants \
        --keep-original-ac \
        -sn ~{sample_ID} \
        -select "QUAL > 20.0 && DP > 5" 
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.lowAC.lowAF.HQ.~{sample_ID}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.lowAC.lowAF.HQ.~{sample_ID}.vcf.gz.tbi"
  }
}

############
### Filter to retrieve High impact variants
############
task HimpactFilter {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String sample_ID
    String base_output_name
  }

command <<<
     java -Xmx8g -jar ~{GATK} \
     SelectVariants \
        -R ~{ref_fasta} \
        -V ~{input_vcf} \
        -O "~{base_output_name}.lowAC.lowAF.HQ.Himpact.~{sample_ID}.vcf.gz" \
        -select "Func.ensGene == 'exonic' || Func.refGene == 'exonic' || dbscSNV_ADA_SCORE > 0.6 || dbscSNV_RF_SCORE > 0.6 || CLINSIG  == 'Pathogenic/Likely_pathogenic' || CLINSIG  == 'Likely_pathogenic' || CLINSIG  == 'Pathogenic'"
##       || dpsi_zscore < -2.0 || dpsi_max_tissue < -2.0
##        -select "ExonicFunc.ensGene != 'synonymous_SNV'"
>>>

#-select "M-CAP_score > 0.025"  

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.lowAC.lowAF.HQ.Himpact.~{sample_ID}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.lowAC.lowAF.HQ.Himpact.~{sample_ID}.vcf.gz.tbi"
  }
}

############
### Convert filtered VCF to Table
############
task VCFToTSV {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String sample_ID
    String base_output_name
  }

command <<<
     java -Xmx8g -jar ~{GATK} \
     VariantsToTable \
        -V ~{input_vcf} \
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AD -F DP -F GQ -F GT \
        -F 1000g2015aug_all -F AAChange.ensGene -F AAChange.refGene -F AC_Orig -F AF_Orig -F AF_afr -F AF_ami -F AF_amr -F AF_asj -F AF_eas -F AF_female -F AF_fin -F AF_male -F AF_nfe -F AF_popmax -F AN_Orig -F AS_FS -F AS_MQ \
        -F AS_QD -F CADD_phred -F Classification -F CLINSIG -F CLNDBN -F CLNDSDBID -F Confidence -F DANN_score -F Disease_description -F DP -F ExcessHet -F ExonicFunc.ensGene -F ExonicFunc.refGene -F FATHMM_pred -F FATHMM_score -F Func.ensGene -F Func.refGene -F GERP++_RS -F Gene.ensGene -F Gene.refGene -F Inheritance \
        -F Kaviar_AF -F M-CAP_pred -F M-CAP_score -F MetaLR_pred -F MetaLR_score -F MetaSVM_pred -F MetaSVM_score -F MOI -F MutationAssessor_pred -F MutationAssessor_score -F MutationTaster_pred -F MutationTaster_score -F NEGATIVE_TRAIN_SITE -F POSITIVE_TRAIN_SITE -F Polyphen2_HDIV_pred -F Polyphen2_HDIV_score \
        -F Polyphen2_HVAR_pred -F Polyphen2_HVAR -F Priority -F QD -F SIFT_pred -F SIFT_score -F max_aaf_all -F controls_AF_popmax -F dbscSNV_ADA_SCORE -F dbscSNV_RF_SCORE -F dpsi_max_tissue -F dpsi_zscore -F non_topmed_AF_popmax -F n.PLP.ClinVar -F n.PLP.LoF.ClinVar -F n.PLP.mis.ClinVar -F non_cancer_AF_popmax -F oe.LoF.upper \
        -F pLI -F phastCons100way_vertebrate -F phyloP100way_vertebrate -F phyloP20way_mammalian \
        -O "~{base_output_name}.lowAC.lowAF.HQ.Himpact.~{sample_ID}.tsv"
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_tsv = "~{base_output_name}.lowAC.lowAF.HQ.Himpact.~{sample_ID}.tsv"
  }
}