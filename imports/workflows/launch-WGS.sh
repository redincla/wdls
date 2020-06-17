#!/bin/bash

########################
#    WGS pipeline v1    #
########################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit, Lausanne, Switzerland.
# Contact: Claire Redin <claire.redin@chuv.ch>

# checks for appropriate input
if [ $# -eq 2 ]; then
 sample_list=$1 #full path to list of samples to process from same batch, one sample / line
 cohort_name=$2 #batch name - should be the same as folder name under which FQ.gz files are stored
else
 echo -e "\n\nLaunching WGS parallelizer script by sample\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-WGS.sh [sample_list] [cohort_name]"
 echo "sample_list: full path to list of samples to process, cohort_name: self-explanatory"
 exit 1
fi

### Set local parameters
export BASEDIR=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/${cohort_name}
export WRKDIR=${BASEDIR}/0620

if ! [ -e $WRKDIR ]; then
	mkdir $WRKDIR
fi

cd ${WRKDIR}

#export cromwell_options="${BASEDIR}/cromwell_options.json"

### Prepare folder structure  -> Adapt full map with relative paths, not absolute paths since FQ files will be relocated
while read sample_ID; do
    if ! [ -e $sample_ID ]; then
	    mkdir $sample_ID
    fi
    if ! [ -e ${sample_ID}/QC ]; then
	    mkdir ${sample_ID}/QC
    fi
    if ! [ -e ${sample_ID}/Raw ]; then
	    mkdir ${sample_ID}/Raw
    fi
    if ! [ -e ${sample_ID}/Processed ]; then
	    mkdir ${sample_ID}/Processed
    fi

    mv ${BASEDIR}/${sample_ID}*fastq.gz ${WRKDIR}/${sample_ID}/Raw/
    mv ${BASEDIR}/${sample_ID}.full_map.tsv ${WRKDIR}/${sample_ID}/Raw/
    sed -i "s+${BASEDIR}+${WRKDIR}/${sample_ID}/Raw+g" ${WRKDIR}/${sample_ID}/Raw/${sample_ID}.full_map.tsv
done < ${sample_list}

# Make template json for inputs
if ! [ -e ${WRKDIR}/json ]; then
	    mkdir ${WRKDIR}/json
fi
# cd /home/credin/scratch/WGS/wdls \
# && womtool inputs WholeGenomeGermlineSingleSample.wdl \
# > ${WRKDIR}/json/template.input.json \
# && cd -

# Write input json for each sample
while read sample_ID; do
  touch ${WRKDIR}/json/${sample_ID}.WGS.input.json
  cat <<EOF > ${WRKDIR}/json/${sample_ID}.WGS.input.json
{
  "WholeGenomeGermlineSingleSample.full_map": "${WRKDIR}/${sample_ID}/Raw/${sample_ID}.full_map.tsv", 
  "WholeGenomeGermlineSingleSample.make_fofn": "true",
  "WholeGenomeGermlineSingleSample.base_file_name": "${sample_ID}",
  "WholeGenomeGermlineSingleSample.final_vcf_base_name": "${sample_ID}_PLUMBING",

  "WholeGenomeGermlineSingleSample.SAMTOOLS": "/software/UHTS/Analysis/samtools/1.10/bin/samtools",
  "WholeGenomeGermlineSingleSample.PICARD": "/software/UHTS/Analysis/picard-tools/2.21.8/bin/picard.jar",
  "WholeGenomeGermlineSingleSample.BWA": "/software/UHTS/Aligner/bwa/0.7.17/bin/bwa",
  "WholeGenomeGermlineSingleSample.GATK": "/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar",
  "WholeGenomeGermlineSingleSample.GATK3": "/software/UHTS/Analysis/GenomeAnalysisTK/3.8.1.0.gf15c1c3ef/bin/GenomeAnalysisTK.jar",
  "WholeGenomeGermlineSingleSample.VerifyBamID": "/software/UHTS/Analysis/VerifyBamID/1.0.6/bin/VerifyBamID",

  "WholeGenomeGermlineSingleSample.ref_fasta": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta",
  "WholeGenomeGermlineSingleSample.ref_index": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.fai",
  "WholeGenomeGermlineSingleSample.ref_dict": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.dict",
  "WholeGenomeGermlineSingleSample.ref_ann": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.64.ann",
  "WholeGenomeGermlineSingleSample.ref_sa": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.64.sa",
  "WholeGenomeGermlineSingleSample.ref_bwt": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.64.bwt",
  "WholeGenomeGermlineSingleSample.ref_pac": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.64.pac",
  "WholeGenomeGermlineSingleSample.ref_amb": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.64.amb",
  "WholeGenomeGermlineSingleSample.ref_alt": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta.64.alt",
  "WholeGenomeGermlineSingleSample.contamination_sites_bed": "/data/chuv/LABO/redin/references/hg38/1000g.phase3.100k.b38.vcf.gz.dat.bed",
  "WholeGenomeGermlineSingleSample.contamination_sites_mu": "/data/chuv/LABO/redin/references/hg38/1000g.phase3.100k.b38.vcf.gz.dat.mu",
  "WholeGenomeGermlineSingleSample.contamination_sites_ud": "/data/chuv/LABO/redin/references/hg38/1000g.phase3.100k.b38.vcf.gz.dat.UD",
  "WholeGenomeGermlineSingleSample.disable_sanity_check": "true",

  "WholeGenomeGermlineSingleSample.wgs_coverage_interval_list": "/data/chuv/LABO/redin/references/hg38/wgs_coverage_regions.hg38.interval_list",
  "WholeGenomeGermlineSingleSample.evaluation_interval_list": "/data/chuv/LABO/redin/references/hg38/wgs_evaluation_regions.hg38.interval_list",
  "WholeGenomeGermlineSingleSample.calling_interval_list": "/data/chuv/LABO/redin/references/hg38/wgs_calling_regions.hg38.interval_list",
  "WholeGenomeGermlineSingleSample.break_bands_at_multiples_of": "100000",
  "WholeGenomeGermlineSingleSample.haplotype_scatter_count": "10",
  "WholeGenomeGermlineSingleSample.make_bamout": "false",

  "WholeGenomeGermlineSingleSample.known_indels_sites_vcfs": [
    "/data/chuv/LABO/redin/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
  ],
  "WholeGenomeGermlineSingleSample.known_indels_sites_indices": [
    "/data/chuv/LABO/redin/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
    "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
  ],
  "WholeGenomeGermlineSingleSample.dbsnp_vcf": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf",
  "WholeGenomeGermlineSingleSample.dbsnp_vcf_index": "/data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
}
EOF
done < ${sample_list}

# Write cromwell options with final output directory for each sample
while read sample_ID; do
  touch ${WRKDIR}/json/${sample_ID}.WGS.options.json
  cat <<EOF > ${WRKDIR}/json/${sample_ID}.WGS.options.json
{
	"final_workflow_outputs_dir": "${WRKDIR}/${sample_ID}/Processed",
  "use_relative_output_paths": true
}
EOF
done < ${sample_list}

### Run first 3 samples (as a test)
# Submit
#while read sample_ID; do
#  sleep 5
#  cd ${WRKDIR}
#  sbatch /home/credin/scratch/WGS/wdls/submit_WGSinglesample.sh \
#    "${WRKDIR}/json/${sample_ID}.WGS.input.json" \
#    "${WRKDIR}/json/${sample_ID}.WGS.options.json"
# done < ${sample_list}

### Run first 3 samples (as a test)
while read sample_ID; do
cd ${WRKDIR}
touch ${WRKDIR}/json/${sample_ID}.WGS.script-submit
echo "#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=WGS
#SBATCH -t 3-00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --output=slurm.%N.%j.log
#module load Utility/cromwell/47

module add Development/java/1.8.0_232
module add R/3.6.1
module add UHTS/Analysis/samtools/1.10
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module add UHTS/Analysis/picard-tools/2.21.8

java -Dconfig.file=/home/credin/.cromwell.conf -jar /software/Utility/cromwell/47/bin/cromwell-47.jar  run /home/credin/scratch/WGS/wdls/WholeGenomeGermlineSingleSample.wdl -i ${WRKDIR}/json/${sample_ID}.WGS.input.json -o ${WRKDIR}/json/${sample_ID}.WGS.options.json" > ${WRKDIR}/json/${sample_ID}.WGS.script-submit

sbatch ${WRKDIR}/json/${sample_ID}.WGS.script-submit

done < ${sample_list}