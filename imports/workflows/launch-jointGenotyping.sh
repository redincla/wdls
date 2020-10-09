#!/bin/bash

####################################
#    JointGenotyping pipeline v1    #
####################################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit, Lausanne, Switzerland.
# Contact: Claire Redin <claire.redin@chuv.ch>

# checks for appropriate input
if [ $# -eq 2 ]; then
 full_map=$1  #col1: sample name, col2: path to gvcf, col3: path to gvcf index
 cohort_name=$2

else
 echo -e "\n\nLaunching jointGenotyping script by sample\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-jointGenotyping.sh [full_map] [cohort_name] "
 echo "full_map: full path to file with 3 tab-delimited columns. Col 1: sample name, col 2: path to gvcf, col 3: path to gvcf index"
 echo "cohort_name: string, name of processed cohort"
 exit 1
fi

### Set local parameters
export BASEDIR=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ
export WRKDIR=${BASEDIR}/${cohort_name}

cut -f 1,2 ${full_map} > $WRKDIR/sample_name_map

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
    touch ${WRKDIR}/json/${cohort_name}.JointGenotyping.input.json
    cat <<EOF > ${WRKDIR}/json/${cohort_name}.JointGenotyping.input.json
{
  "JointGenotyping.full_map": "${full_map}",
  "JointGenotyping.sample_name_map": "$WRKDIR/sample_name_map",
  "JointGenotyping.sample_num_threshold": "50",

  "JointGenotyping.GATK": "/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar",
  "JointGenotyping.PICARD": "/software/UHTS/Analysis/picard-tools/2.21.8/bin/picard.jar",
  "JointGenotyping.tabix": "/software/UHTS/Analysis/EPACTS/3.2.6/bin/tabix",

  "JointGenotyping.unpadded_intervals_file": "/data/chuv/LABO/redin/references/hg38/hg38.even.handcurated.20k.intervals",
  "JointGenotyping.ref_fasta": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta",
  "JointGenotyping.ref_index": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.fai",
  "JointGenotyping.ref_dict": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.dict",

  "JointGenotyping.workspace_dir_name": "genomicsdb",
  "JointGenotyping.callset_name": "${cohort_name}",

  "JointGenotyping.top_level_scatter_count": "40",

  "JointGenotyping.dbsnp_vcf": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf",
  "JointGenotyping.dbsnp_vcf_index": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.idx",

  "JointGenotyping.indel_recalibration_tranche_values": ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"],
  "JointGenotyping.indel_recalibration_annotation_values": ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"],
  "JointGenotyping.snp_recalibration_tranche_values": ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ],
  "JointGenotyping.snp_recalibration_annotation_values": ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"],

  "JointGenotyping.mills_resource_vcf": "/data/chuv/LABO/redin/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
  "JointGenotyping.mills_resource_vcf_index": "/data/chuv/LABO/redin/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
  "JointGenotyping.axiomPoly_resource_vcf": "/data/chuv/LABO/redin/references/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
  "JointGenotyping.axiomPoly_resource_vcf_index": "/data/chuv/LABO/redin/references/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
  "JointGenotyping.hapmap_resource_vcf": "/data/chuv/LABO/redin/references/hg38/hapmap_3.3.hg38.vcf.gz",
  "JointGenotyping.hapmap_resource_vcf_index": "/data/chuv/LABO/redin/references/hg38/hapmap_3.3.hg38.vcf.gz.tbi",
  "JointGenotyping.omni_resource_vcf": "/data/chuv/LABO/redin/references/hg38/1000G_omni2.5.hg38.vcf.gz",
  "JointGenotyping.omni_resource_vcf_index": "/data/chuv/LABO/redin/references/hg38/1000G_omni2.5.hg38.vcf.gz.tbi",
  "JointGenotyping.one_thousand_genomes_resource_vcf": "/data/chuv/LABO/redin/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
  "JointGenotyping.one_thousand_genomes_resource_vcf_index": "/data/chuv/LABO/redin/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
 
  "JointGenotyping.snp_filter_level": "99.7",
  "JointGenotyping.indel_filter_level": "99.0",

  "JointGenotyping.eval_interval_list": "/data/chuv/LABO/redin/references/hg38/wgs_evaluation_regions.hg38.interval_list",

  "JointGenotyping.haplotype_database": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.haplotype_database.txt",

  "JointGenotyping.use_allele_specific_annotations": "false",
  "JointGenotyping.cross_check_fingerprints": "true",
  "JointGenotyping.excess_het_threshold": "54.69"
}
EOF

# Write cromwell options with final output directory
  touch ${WRKDIR}/json/${cohort_name}.JointGenotyping.options.json
  cat <<EOF > ${WRKDIR}/json/${cohort_name}.JointGenotyping.options.json
{
    "final_workflow_outputs_dir": "$WRKDIR/JointCalling",
    "use_relative_output_paths": true
}
EOF

### Launch jointGenotyping script
cd ${WRKDIR}
touch ${WRKDIR}/json/${cohort_name}.JointGenotyping.script-submit
echo "#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=jointCalling.${cohort_name}
#SBATCH -t 6-00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --output=${WRKDIR}/json/slurm.${cohort_name}.%N.%j.log

module add Development/java/1.8.0_232
module add R/3.6.1
module add UHTS/Analysis/samtools/1.10
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module add UHTS/Analysis/picard-tools/2.21.8

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/imports/workflows/JointGenotyping.wdl -i ${WRKDIR}/json/${cohort_name}.JointGenotyping.input.json -o ${WRKDIR}/json/${cohort_name}.JointGenotyping.options.json" > ${WRKDIR}/json/${cohort_name}.JointGenotyping.script-submit

sbatch ${WRKDIR}/json/${cohort_name}.JointGenotyping.script-submit