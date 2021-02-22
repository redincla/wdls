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
    call SplitbySample { 
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
            input_vcf = SplitbySample.output_vcf,
            input_vcf_index = SplitbySample.output_vcf_index,
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
    gzip -cd ~{input_vcf} | sed 's@dbscSNV_ADA_SCORE,Number=.,Type=String@dbscSNV_ADA_SCORE,Number=1,Type=Float@g' | sed 's@dbscSNV_RF_SCORE,Number=.,Type=String@dbscSNV_RF_SCORE,Number=1,Type=Float@g' | sed 's@dpsi_max_tissue,Number=.,Type=String@dpsi_max_tissue,Number=1,Type=Float@g' | sed 's@dpsi_zscore,Number=.,Type=String@dpsi_zscore,Number=1,Type=Float@g' | sed 's@;dbscSNV_ADA_SCORE=\.;@;dbscSNV_ADA_SCORE=NaN;@g'  | sed 's@;dbscSNV_RF_SCORE=\.;@;dbscSNV_RF_SCORE=NaN;@g' | sed 's@;dpsi_max_tissue=\.;@;dpsi_max_tissue=NaN;@g' | sed 's@;dpsi_zscore=\.;@;dpsi_zscore=NaN;@g' | sed 's@;max_aaf_all=\.;@;max_aaf_all=NaN;@g' | bgzip -c > tmp.vcf.gz
    tabix tmp.vcf.gz
    java -Xmx8g -jar ~{GATK} \
    SelectVariants \
    -R ~{ref_fasta} \
    -V tmp.vcf.gz \
    -O ~{base_output_name}.lowAC.lowAF.vcf.gz \
    -select "max_aaf_all < 0.01 || max_aaf_all == 'NaN'"
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
### Split by sample && filter on GQ/DP -> already flagged in GenotypeRefinement workflow. Redundant.
############
task SplitbySample {
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
        -O "~{base_output_name}.lowAC.lowAF.~{sample_ID}.vcf.gz" \
        --exclude-non-variants \
        --keep-original-ac \
        -sn ~{sample_ID}
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.lowAC.lowAF.~{sample_ID}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.lowAC.lowAF.~{sample_ID}.vcf.gz.tbi"
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
        -O "~{base_output_name}.lowAC.lowAF.Himpact.~{sample_ID}.vcf.gz" \
        -select "Func.ensGene == 'exonic' || Func.refGene == 'exonic' || dbscSNV_ADA_SCORE > 0.6 || dbscSNV_RF_SCORE > 0.6 || ClinVar_Sign  == 'Pathogenic/Likely_pathogenic' || ClinVar_Sign  == 'Likely_pathogenic' || ClinVar_Sign  == 'Pathogenic'"
##       || dpsi_zscore < -2.0 || dpsi_max_tissue < -2.0
##        -select "ExonicFunc.ensGene != 'synonymous_SNV'"
>>>

#-select "M-CAP_score > 0.025"  

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_vcf = "~{base_output_name}.lowAC.lowAF.Himpact.~{sample_ID}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.lowAC.lowAF.Himpact.~{sample_ID}.vcf.gz.tbi"
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
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -GF AD -F DP -GF GQ -GF GT -F AC_Orig -F AF_Orig -F AN_Orig -F AS_FS -F AS_MQ -F AS_QD -F QD -F DP -F ExcessHet -F NEGATIVE_TRAIN_SITE -F POSITIVE_TRAIN_SITE -F SB \
        -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F NSAMPLES -F NCALLED \
        -F Gene.ensGene -F Gene.refGene -F Func.ensGene -F Func.refGene -F ExonicFunc.ensGene -F ExonicFunc.refGene -F AAChange.ensGene -F AAChange.refGene \
        -F gnomadv3_AF -F gnomadv3_AF_raw -F gnomadv3_AF_male -F gnomadv3_AF_female -F gnomadv3_AF_afr -F gnomadv3_AF_ami -F gnomadv3_AF_amr -F gnomadv3_AF_asj -F gnomadv3_AF_eas -F gnomadv3_AF_fin -F gnomadv3_AF_nfe -F gnomadv3_AF_oth -F gnomadv3_AF_sas -F Kaviar_AF -F 1000g2015aug_all -F ClinVar_AF_EXAC -F ClinVar_AF_TGP -F max_aaf_all \
        -F ClinVar_Sign -F CLINSIG -F ClinVar_Sign_Conflict -F ClinVar_RevStatus -F ClinVar_Disease_name -F ClinGen_Disease_name -F ClinGen_MOI -F ClinGen_Classification -F MIM_Disease -F MIM_Disease_name -F MIM_gene -F DDD_mutation_consequence -F DDD_MOI -F DDD_Classification -F ACMG -F MEDP_MOI -F MEDP_Confidence -F n.PLP.ClinVar -F n.PLP.LoF.ClinVar -F n.PLP.mis.ClinVar -F IFNpathway_gene_name -F IFNpathway_gene_id \
        -F DamagePredCount -F SIFT_pred -F SIFT4G_pred -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F LRT_pred -F MutationTaster_pred -F MutationAssessor_pred -F FATHMM_pred -F PROVEAN_pred -F VEST4_score -F MetaSVM_pred -F MetaLR_pred -F M-CAP_pred -F REVEL_score -F MutPred_score -F MVP_score -F MPC_score \
        -F DEOGEN2_pred -F BayesDel_addAF_pred -F BayesDel_noAF_pred -F ClinPred_pred -F LIST-S2_pred -F CADD_raw -F CADD_phred -F DANN_score -F fathmm-MKL_coding_pred -F fathmm-XF_coding_pred -F Eigen-raw_coding -F Eigen-phred_coding -F LOEUF_bin -F GnomAD_pLI -F ExAC_pLI \
        -F PrimateAI_pred -F Eigen-PC-raw_coding -F Eigen-PC-phred_coding -F GenoCanyon_score -F LINSIGHT -F GERP++_NR -F GERP++_RS -F phyloP100way_vertebrate -F phyloP30way_mammalian -F phyloP17way_primate -F phastCons100way_vertebrate -F phastCons30way_mammalian -F phastCons17way_primate -F oe.LoF.upper -F pLI \
        -F bStatistic -F Interpro_domain -F GTEx_V8_gene -F GTEx_V8_tissue -F integrated_fitCons_score \
        -F dbscSNV_ADA_SCORE -F dbscSNV_RF_SCORE -F dpsi_max_tissue -F dpsi_zscore  \
        -O "~{base_output_name}.lowAC.lowAF.Himpact.~{sample_ID}.tsv"
>>>

  runtime {
  cpus: "1"
	requested_memory_mb_per_core: "9000" 
  }

  output {
    File output_tsv = "~{base_output_name}.lowAC.lowAF.Himpact.~{sample_ID}.tsv"
  }
}