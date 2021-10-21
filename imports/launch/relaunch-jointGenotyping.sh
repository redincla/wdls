#!/bin/bash

####################################
#    JointGenotyping pipeline v1    #
####################################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit, Lausanne, Switzerland.
# Contact: Claire Redin <claire.redin@chuv.ch>

# checks for appropriate input
if [ $# -eq 3 ]; then
 full_map=$1  #col1: sample name, col2: path to gvcf, col3: path to gvcf index
 cohort_name=$2
 WRKDIR=$3
else
 echo -e "\n\nLaunching jointGenotyping script by sample\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-jointGenotyping.sh [full_map] [cohort_name] [output directory]"
 echo "full_map: full path to file with 3 tab-delimited columns. Col 1: sample name, col 2: path to gvcf, col 3: path to gvcf index"
 echo "cohort_name: string, name of processed cohort"
 echo "output directory: full path"
 exit 1
fi

### Set local parameters
export WRKDIR

cut -f 1,2 ${full_map} > $WRKDIR/sample_name_map
sample_number=$(wc -l ${full_map} | awk '{print $1}')

if ! [ -e $WRKDIR ]; then
	mkdir $WRKDIR
fi

cd ${WRKDIR}

### Prepare folder structure
if ! [ -e $WRKDIR/JointCalling ]; then
	mkdir $WRKDIR/JointCalling
fi

if ! [ -e ${WRKDIR}/json ]; then
	    mkdir ${WRKDIR}/json
fi

# Write input json
    touch ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.input.json

    shard_number=$( echo $sample_number / 3 | bc )  ### calculate number of shards needed based on total # of input genomes
    cat <<EOF > ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.input.json
{
  "JointGenotyping.full_map": "${full_map}",
  "JointGenotyping.sample_name_map": "$WRKDIR/sample_name_map",
  "JointGenotyping.sample_num_threshold": "50",

  "JointGenotyping.GATK": "/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar",
  "JointGenotyping.PICARD": "/software/UHTS/Analysis/picard-tools/2.21.8/bin/picard.jar",
  "JointGenotyping.tabix": "/software/UHTS/Analysis/EPACTS/3.2.6/bin/tabix",

  "JointGenotyping.unpadded_intervals_file": "/home/credin/refs/references/hg38/hg38.even.handcurated.20k.intervals",
  "JointGenotyping.ref_fasta": "/home/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta",
  "JointGenotyping.ref_index": "/home/credin/refs/references/hg38/Homo_sapiens_assembly38.fasta.fai",
  "JointGenotyping.ref_dict": "/home/credin/refs/references/hg38/Homo_sapiens_assembly38.dict",

  "JointGenotyping.workspace_dir_name": "genomicsdb",
  "JointGenotyping.callset_name": "${cohort_name}",

  "JointGenotyping.top_level_scatter_count": "${shard_number}",

  "JointGenotyping.dbsnp_vcf": "/home/credin/refs/references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf",
  "JointGenotyping.dbsnp_vcf_index": "/home/credin/refs/references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.idx",

  "JointGenotyping.indel_recalibration_tranche_values": ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"],
  "JointGenotyping.indel_recalibration_annotation_values": ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"],
  "JointGenotyping.snp_recalibration_tranche_values": ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ],
  "JointGenotyping.snp_recalibration_annotation_values": ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"],

  "JointGenotyping.mills_resource_vcf": "/home/credin/refs/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
  "JointGenotyping.mills_resource_vcf_index": "/home/credin/refs/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
  "JointGenotyping.axiomPoly_resource_vcf": "/home/credin/refs/references/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
  "JointGenotyping.axiomPoly_resource_vcf_index": "/home/credin/refs/references/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
  "JointGenotyping.hapmap_resource_vcf": "/home/credin/refs/references/hg38/hapmap_3.3.hg38.vcf.gz",
  "JointGenotyping.hapmap_resource_vcf_index": "/home/credin/refs/references/hg38/hapmap_3.3.hg38.vcf.gz.tbi",
  "JointGenotyping.omni_resource_vcf": "/home/credin/refs/references/hg38/1000G_omni2.5.hg38.vcf.gz",
  "JointGenotyping.omni_resource_vcf_index": "/home/credin/refs/references/hg38/1000G_omni2.5.hg38.vcf.gz.tbi",
  "JointGenotyping.one_thousand_genomes_resource_vcf": "/home/credin/refs/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
  "JointGenotyping.one_thousand_genomes_resource_vcf_index": "/home/credin/refs/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
  "JointGenotyping.one_thousand_genomes_resource_NO_MULTIALLELIC_vcf": "/home/credin/refs/references/hg38/1000G_phase1.snps.high_confidence.NO_MULTIALLELIC.hg38.vcf.gz",
  "JointGenotyping.one_thousand_genomes_resource_NO_MULTIALLELIC_vcf_index": "/home/credin/refs/references/hg38/1000G_phase1.snps.high_confidence.NO_MULTIALLELIC.hg38.vcf.gz.tbi",

  "JointGenotyping.snp_filter_level": "99.7",
  "JointGenotyping.indel_filter_level": "99.0",

  "JointGenotyping.eval_interval_list": "/home/credin/refs/references/hg38/wgs_evaluation_regions.hg38.interval_list",

  "JointGenotyping.haplotype_database": "/home/credin/refs/references/hg38/Homo_sapiens_assembly38.haplotype_database.txt",

  "JointGenotyping.use_allele_specific_annotations": "false",
  "JointGenotyping.cross_check_fingerprints": "true",
  "JointGenotyping.excess_het_threshold": "54.69",

  "JointGenotyping.genotyped_vcf": [
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-61/execution/210708_TACG_CoPRAC_HEV_VZV.61.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-10/execution/210708_TACG_CoPRAC_HEV_VZV.10.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-34/execution/210708_TACG_CoPRAC_HEV_VZV.34.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-60/execution/210708_TACG_CoPRAC_HEV_VZV.60.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-45/execution/210708_TACG_CoPRAC_HEV_VZV.45.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-66/execution/210708_TACG_CoPRAC_HEV_VZV.66.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-2/execution/210708_TACG_CoPRAC_HEV_VZV.2.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-26/execution/210708_TACG_CoPRAC_HEV_VZV.26.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-36/execution/210708_TACG_CoPRAC_HEV_VZV.36.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-50/execution/210708_TACG_CoPRAC_HEV_VZV.50.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-55/execution/210708_TACG_CoPRAC_HEV_VZV.55.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-56/execution/210708_TACG_CoPRAC_HEV_VZV.56.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-58/execution/210708_TACG_CoPRAC_HEV_VZV.58.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-53/execution/210708_TACG_CoPRAC_HEV_VZV.53.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-20/execution/210708_TACG_CoPRAC_HEV_VZV.20.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-39/execution/210708_TACG_CoPRAC_HEV_VZV.39.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-22/execution/210708_TACG_CoPRAC_HEV_VZV.22.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-59/execution/210708_TACG_CoPRAC_HEV_VZV.59.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-51/execution/210708_TACG_CoPRAC_HEV_VZV.51.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-49/execution/210708_TACG_CoPRAC_HEV_VZV.49.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-37/execution/210708_TACG_CoPRAC_HEV_VZV.37.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-44/execution/210708_TACG_CoPRAC_HEV_VZV.44.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-19/execution/210708_TACG_CoPRAC_HEV_VZV.19.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-41/execution/210708_TACG_CoPRAC_HEV_VZV.41.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-18/execution/210708_TACG_CoPRAC_HEV_VZV.18.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-5/execution/210708_TACG_CoPRAC_HEV_VZV.5.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-11/execution/210708_TACG_CoPRAC_HEV_VZV.11.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-8/execution/210708_TACG_CoPRAC_HEV_VZV.8.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-57/execution/210708_TACG_CoPRAC_HEV_VZV.57.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-32/execution/210708_TACG_CoPRAC_HEV_VZV.32.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-48/execution/210708_TACG_CoPRAC_HEV_VZV.48.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-14/execution/210708_TACG_CoPRAC_HEV_VZV.14.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-17/execution/210708_TACG_CoPRAC_HEV_VZV.17.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-46/execution/210708_TACG_CoPRAC_HEV_VZV.46.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-16/execution/210708_TACG_CoPRAC_HEV_VZV.16.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-6/execution/210708_TACG_CoPRAC_HEV_VZV.6.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-13/execution/210708_TACG_CoPRAC_HEV_VZV.13.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-23/execution/210708_TACG_CoPRAC_HEV_VZV.23.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-30/execution/210708_TACG_CoPRAC_HEV_VZV.30.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-28/execution/210708_TACG_CoPRAC_HEV_VZV.28.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-52/execution/210708_TACG_CoPRAC_HEV_VZV.52.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-47/execution/210708_TACG_CoPRAC_HEV_VZV.47.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-43/execution/210708_TACG_CoPRAC_HEV_VZV.43.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-1/execution/210708_TACG_CoPRAC_HEV_VZV.1.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-40/execution/210708_TACG_CoPRAC_HEV_VZV.40.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-15/execution/210708_TACG_CoPRAC_HEV_VZV.15.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-62/execution/210708_TACG_CoPRAC_HEV_VZV.62.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-7/execution/210708_TACG_CoPRAC_HEV_VZV.7.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-33/execution/210708_TACG_CoPRAC_HEV_VZV.33.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-65/execution/210708_TACG_CoPRAC_HEV_VZV.65.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-4/execution/210708_TACG_CoPRAC_HEV_VZV.4.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-38/execution/210708_TACG_CoPRAC_HEV_VZV.38.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-63/execution/210708_TACG_CoPRAC_HEV_VZV.63.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-42/execution/210708_TACG_CoPRAC_HEV_VZV.42.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-12/execution/210708_TACG_CoPRAC_HEV_VZV.12.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-9/execution/210708_TACG_CoPRAC_HEV_VZV.9.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-0/execution/210708_TACG_CoPRAC_HEV_VZV.0.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-35/execution/210708_TACG_CoPRAC_HEV_VZV.35.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-3/execution/210708_TACG_CoPRAC_HEV_VZV.3.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-54/execution/210708_TACG_CoPRAC_HEV_VZV.54.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-31/execution/210708_TACG_CoPRAC_HEV_VZV.31.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-27/execution/210708_TACG_CoPRAC_HEV_VZV.27.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-25/execution/210708_TACG_CoPRAC_HEV_VZV.25.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-24/execution/210708_TACG_CoPRAC_HEV_VZV.24.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-21/execution/210708_TACG_CoPRAC_HEV_VZV.21.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-29/execution/210708_TACG_CoPRAC_HEV_VZV.29.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-64/execution/210708_TACG_CoPRAC_HEV_VZV.64.vcf.gz"],

"JointGenotyping.genotyped_vcf_index": [
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-61/execution/210708_TACG_CoPRAC_HEV_VZV.61.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-10/execution/210708_TACG_CoPRAC_HEV_VZV.10.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-34/execution/210708_TACG_CoPRAC_HEV_VZV.34.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-60/execution/210708_TACG_CoPRAC_HEV_VZV.60.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-45/execution/210708_TACG_CoPRAC_HEV_VZV.45.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-66/execution/210708_TACG_CoPRAC_HEV_VZV.66.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-2/execution/210708_TACG_CoPRAC_HEV_VZV.2.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-26/execution/210708_TACG_CoPRAC_HEV_VZV.26.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-36/execution/210708_TACG_CoPRAC_HEV_VZV.36.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-50/execution/210708_TACG_CoPRAC_HEV_VZV.50.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-55/execution/210708_TACG_CoPRAC_HEV_VZV.55.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-56/execution/210708_TACG_CoPRAC_HEV_VZV.56.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-58/execution/210708_TACG_CoPRAC_HEV_VZV.58.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-53/execution/210708_TACG_CoPRAC_HEV_VZV.53.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-20/execution/210708_TACG_CoPRAC_HEV_VZV.20.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-39/execution/210708_TACG_CoPRAC_HEV_VZV.39.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-22/execution/210708_TACG_CoPRAC_HEV_VZV.22.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-59/execution/210708_TACG_CoPRAC_HEV_VZV.59.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-51/execution/210708_TACG_CoPRAC_HEV_VZV.51.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-49/execution/210708_TACG_CoPRAC_HEV_VZV.49.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-37/execution/210708_TACG_CoPRAC_HEV_VZV.37.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-44/execution/210708_TACG_CoPRAC_HEV_VZV.44.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-19/execution/210708_TACG_CoPRAC_HEV_VZV.19.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-41/execution/210708_TACG_CoPRAC_HEV_VZV.41.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-18/execution/210708_TACG_CoPRAC_HEV_VZV.18.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-5/execution/210708_TACG_CoPRAC_HEV_VZV.5.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-11/execution/210708_TACG_CoPRAC_HEV_VZV.11.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-8/execution/210708_TACG_CoPRAC_HEV_VZV.8.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-57/execution/210708_TACG_CoPRAC_HEV_VZV.57.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-32/execution/210708_TACG_CoPRAC_HEV_VZV.32.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-48/execution/210708_TACG_CoPRAC_HEV_VZV.48.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-14/execution/210708_TACG_CoPRAC_HEV_VZV.14.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-17/execution/210708_TACG_CoPRAC_HEV_VZV.17.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-46/execution/210708_TACG_CoPRAC_HEV_VZV.46.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-16/execution/210708_TACG_CoPRAC_HEV_VZV.16.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-6/execution/210708_TACG_CoPRAC_HEV_VZV.6.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-13/execution/210708_TACG_CoPRAC_HEV_VZV.13.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-23/execution/210708_TACG_CoPRAC_HEV_VZV.23.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-30/execution/210708_TACG_CoPRAC_HEV_VZV.30.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-28/execution/210708_TACG_CoPRAC_HEV_VZV.28.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-52/execution/210708_TACG_CoPRAC_HEV_VZV.52.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-47/execution/210708_TACG_CoPRAC_HEV_VZV.47.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-43/execution/210708_TACG_CoPRAC_HEV_VZV.43.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-1/execution/210708_TACG_CoPRAC_HEV_VZV.1.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-40/execution/210708_TACG_CoPRAC_HEV_VZV.40.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-15/execution/210708_TACG_CoPRAC_HEV_VZV.15.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-62/execution/210708_TACG_CoPRAC_HEV_VZV.62.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-7/execution/210708_TACG_CoPRAC_HEV_VZV.7.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-33/execution/210708_TACG_CoPRAC_HEV_VZV.33.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-65/execution/210708_TACG_CoPRAC_HEV_VZV.65.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-4/execution/210708_TACG_CoPRAC_HEV_VZV.4.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-38/execution/210708_TACG_CoPRAC_HEV_VZV.38.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-63/execution/210708_TACG_CoPRAC_HEV_VZV.63.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-42/execution/210708_TACG_CoPRAC_HEV_VZV.42.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-12/execution/210708_TACG_CoPRAC_HEV_VZV.12.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-9/execution/210708_TACG_CoPRAC_HEV_VZV.9.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-0/execution/210708_TACG_CoPRAC_HEV_VZV.0.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-35/execution/210708_TACG_CoPRAC_HEV_VZV.35.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-3/execution/210708_TACG_CoPRAC_HEV_VZV.3.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-54/execution/210708_TACG_CoPRAC_HEV_VZV.54.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-31/execution/210708_TACG_CoPRAC_HEV_VZV.31.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-27/execution/210708_TACG_CoPRAC_HEV_VZV.27.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-25/execution/210708_TACG_CoPRAC_HEV_VZV.25.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-24/execution/210708_TACG_CoPRAC_HEV_VZV.24.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-21/execution/210708_TACG_CoPRAC_HEV_VZV.21.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-29/execution/210708_TACG_CoPRAC_HEV_VZV.29.vcf.gz.tbi",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-GenotypeGVCFs/shard-64/execution/210708_TACG_CoPRAC_HEV_VZV.64.vcf.gz.tbi"],

"JointGenotyping.sites_only_vcfs": [
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-61/execution/210708_TACG_CoPRAC_HEV_VZV.61.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-66/execution/210708_TACG_CoPRAC_HEV_VZV.66.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-45/execution/210708_TACG_CoPRAC_HEV_VZV.45.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-60/execution/210708_TACG_CoPRAC_HEV_VZV.60.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-34/execution/210708_TACG_CoPRAC_HEV_VZV.34.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-37/execution/210708_TACG_CoPRAC_HEV_VZV.37.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-50/execution/210708_TACG_CoPRAC_HEV_VZV.50.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-19/execution/210708_TACG_CoPRAC_HEV_VZV.19.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-47/execution/210708_TACG_CoPRAC_HEV_VZV.47.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-44/execution/210708_TACG_CoPRAC_HEV_VZV.44.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-32/execution/210708_TACG_CoPRAC_HEV_VZV.32.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-57/execution/210708_TACG_CoPRAC_HEV_VZV.57.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-2/execution/210708_TACG_CoPRAC_HEV_VZV.2.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-13/execution/210708_TACG_CoPRAC_HEV_VZV.13.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-6/execution/210708_TACG_CoPRAC_HEV_VZV.6.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-10/execution/210708_TACG_CoPRAC_HEV_VZV.10.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-53/execution/210708_TACG_CoPRAC_HEV_VZV.53.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-33/execution/210708_TACG_CoPRAC_HEV_VZV.33.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-22/execution/210708_TACG_CoPRAC_HEV_VZV.22.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-20/execution/210708_TACG_CoPRAC_HEV_VZV.20.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-59/execution/210708_TACG_CoPRAC_HEV_VZV.59.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-39/execution/210708_TACG_CoPRAC_HEV_VZV.39.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-49/execution/210708_TACG_CoPRAC_HEV_VZV.49.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-3/execution/210708_TACG_CoPRAC_HEV_VZV.3.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-5/execution/210708_TACG_CoPRAC_HEV_VZV.5.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-52/execution/210708_TACG_CoPRAC_HEV_VZV.52.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-11/execution/210708_TACG_CoPRAC_HEV_VZV.11.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-51/execution/210708_TACG_CoPRAC_HEV_VZV.51.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-18/execution/210708_TACG_CoPRAC_HEV_VZV.18.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-56/execution/210708_TACG_CoPRAC_HEV_VZV.56.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-41/execution/210708_TACG_CoPRAC_HEV_VZV.41.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-55/execution/210708_TACG_CoPRAC_HEV_VZV.55.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-30/execution/210708_TACG_CoPRAC_HEV_VZV.30.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-26/execution/210708_TACG_CoPRAC_HEV_VZV.26.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-24/execution/210708_TACG_CoPRAC_HEV_VZV.24.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-17/execution/210708_TACG_CoPRAC_HEV_VZV.17.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-28/execution/210708_TACG_CoPRAC_HEV_VZV.28.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-16/execution/210708_TACG_CoPRAC_HEV_VZV.16.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-43/execution/210708_TACG_CoPRAC_HEV_VZV.43.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-23/execution/210708_TACG_CoPRAC_HEV_VZV.23.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-36/execution/210708_TACG_CoPRAC_HEV_VZV.36.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-21/execution/210708_TACG_CoPRAC_HEV_VZV.21.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-42/execution/210708_TACG_CoPRAC_HEV_VZV.42.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-35/execution/210708_TACG_CoPRAC_HEV_VZV.35.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-4/execution/210708_TACG_CoPRAC_HEV_VZV.4.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-58/execution/210708_TACG_CoPRAC_HEV_VZV.58.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-8/execution/210708_TACG_CoPRAC_HEV_VZV.8.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-0/execution/210708_TACG_CoPRAC_HEV_VZV.0.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-7/execution/210708_TACG_CoPRAC_HEV_VZV.7.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-12/execution/210708_TACG_CoPRAC_HEV_VZV.12.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-29/execution/210708_TACG_CoPRAC_HEV_VZV.29.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-31/execution/210708_TACG_CoPRAC_HEV_VZV.31.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-40/execution/210708_TACG_CoPRAC_HEV_VZV.40.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-54/execution/210708_TACG_CoPRAC_HEV_VZV.54.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-9/execution/210708_TACG_CoPRAC_HEV_VZV.9.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-48/execution/210708_TACG_CoPRAC_HEV_VZV.48.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-62/execution/210708_TACG_CoPRAC_HEV_VZV.62.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-14/execution/210708_TACG_CoPRAC_HEV_VZV.14.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-65/execution/210708_TACG_CoPRAC_HEV_VZV.65.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-27/execution/210708_TACG_CoPRAC_HEV_VZV.27.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-1/execution/210708_TACG_CoPRAC_HEV_VZV.1.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-25/execution/210708_TACG_CoPRAC_HEV_VZV.25.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-15/execution/210708_TACG_CoPRAC_HEV_VZV.15.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-46/execution/210708_TACG_CoPRAC_HEV_VZV.46.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-38/execution/210708_TACG_CoPRAC_HEV_VZV.38.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-63/execution/210708_TACG_CoPRAC_HEV_VZV.63.sites_only.vcf.gz",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-MakeSitesOnlyVcf/shard-64/execution/210708_TACG_CoPRAC_HEV_VZV.64.sites_only.vcf.gz"],

"JointGenotyping.unpadded_intervals": [
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0000-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0001-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0002-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0003-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0004-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0005-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0006-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0007-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0008-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0009-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0010-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0011-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0012-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0013-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0014-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0015-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0016-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0017-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0018-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0019-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0020-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0021-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0022-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0023-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0024-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0025-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0026-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0027-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0028-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0029-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0030-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0031-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0032-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0033-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0034-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0035-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0036-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0037-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0038-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0039-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0040-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0041-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0042-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0043-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0044-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0045-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0046-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0047-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0048-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0049-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0050-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0051-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0052-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0053-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0054-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0055-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0056-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0057-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0058-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0059-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0060-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0061-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0062-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0063-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0064-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0065-scattered.interval_list",
"/users/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG/0920/cromwell-executions/JointGenotyping/4b20c928-173f-463e-afa4-6f0440b7107e/call-SplitIntervalList/execution/glob-d928cd0f5fb17b6bd5e635f48c18ccfb/0066-scattered.interval_list"]
}
EOF

# Write cromwell options with final output directory
  touch ${WRKDIR}/json/${cohort_name}.JointGenotyping.options.json
  cat <<EOF > ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.options.json
{
    "final_workflow_outputs_dir": "$WRKDIR/JointCalling",
    "use_relative_output_paths": true
}
EOF

### Launch jointGenotyping script
cd ${WRKDIR}
touch ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.script-submit
echo "#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=jointCalling_relaunch.${cohort_name}
#SBATCH -t 7-00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --output=${WRKDIR}/json/jointCalling_relaunch.slurm.${cohort_name}.%N.%j.log

module add Development/java/1.8.0_232
module add R/3.6.1
module add UHTS/Analysis/samtools/1.10
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module add UHTS/Analysis/picard-tools/2.21.8

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/workflows/JointGenotyping_relaunch.wdl -i ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.input.json -o ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.options.json" > ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.script-submit

sbatch ${WRKDIR}/json/${cohort_name}.JointGenotyping_relaunch.script-submit