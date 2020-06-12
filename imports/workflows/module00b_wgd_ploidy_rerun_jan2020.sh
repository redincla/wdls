#!/usr/bin/env bash

######################
#    gnomAD-SV v3    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Rerun module 00b for corrected WGD & ploidy estimates on gnomAD v3

# Run on local MBP & cromwell via cromshell


### Set local parameters
export BASEDIR=/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD_v3
export WRKDIR=${BASEDIR}/cromshell_local/hg38_module00b_wgd_ploidy_rerun_Jan2020
if ! [ -e $WRKDIR ]; then
	mkdir $WRKDIR
fi
cd ${WRKDIR}
export GATKSV=/Users/rlc/Desktop/Collins/Talkowski/code/gatk-sv-v1
export cromwell_options="${BASEDIR}/gnomad_sv_v3/cromshell/rlc_gnomad_cromshell_options.json"
export cromwell_options_nocache="${BASEDIR}/gnomad_sv_v3/cromshell/rlc_gnomad_cromshell_options.no_call_caching.json"
# export candidate_bins="gs://talkowski-sv-gnomad-resources/rlc_refs/gnomAD-SV_v3.WGD_hg38.candidate_bins.bed.gz"


### Reset directory tree
if ! [ -e ${WRKDIR} ]; then
  mkdir ${WRKDIR}
fi
for SUBDIR in json logs; do
  if ! [ -e ${WRKDIR}/${SUBDIR} ]; then
    mkdir ${WRKDIR}/${SUBDIR}
  fi
done
for SUBDIR in data/mod00b_wgd_ploidy; do
  if ! [ -e ${BASEDIR}/${SUBDIR} ]; then
    mkdir ${BASEDIR}/${SUBDIR}
  fi
done


### Rebuild docker corresponding to WGD/ploidy patch branch
# Checkout fixed branch
cd $GATKSV \
&& git checkout wgd_sample_id_hotfix \
&& git pull
# Use supplied python script to rebuild images
cd $GATKSV/scripts/docker \
&& python build_docker_lib.py \
     --targets sv-base-mini sv-pipeline sv-pipeline-qc \
     --image-tag wgd_sample_id_hotfix-d12bcf \
     --dockerhub-root rlcollins


### Create batches of samples based on master metadata manifest 
# Reset directory, if needed
if [ -e ${WRKDIR}/batches ]; then
  rm -rf ${WRKDIR}/batches
fi
mkdir ${WRKDIR}/batches
# Get sample IDs and coverage file paths
counts_idx=$( zcat ${BASEDIR}/sample_metadata/gnomAD-SV_v3.master_sample_metadata.tsv.gz \
              | head -n1 | sed 's/\t/\n/g' \
              | awk '{ if ($1=="counts_file") print NR }' )
zcat ${BASEDIR}/sample_metadata/gnomAD-SV_v3.master_sample_metadata.tsv.gz \
| cut -f1,${counts_idx} \
> ${WRKDIR}/all_samples_w_cov_paths.tsv
# Clean divisions of 500 samples per batch
~/Desktop/Collins/Talkowski/code/GenomicsToolbox/evenSplitter.R \
  -L 500 \
  --shuffle \
  ${WRKDIR}/all_samples_w_cov_paths.tsv \
  ${WRKDIR}/batches/gnomAD-SV_v3.wgd_ploidy_batch.
n_splits=$( find ${WRKDIR}/batches/ -name 'gnomAD-SV_v3.wgd_ploidy_batch.[0-9]*' | wc -l )
for i in $( seq 1 ${n_splits} ); do
  mv ${WRKDIR}/batches/gnomAD-SV_v3.wgd_ploidy_batch.${i} \
     ${WRKDIR}/batches/gnomAD-SV_v3.wgd_ploidy_batch.batch_${i}.samples_wCovPath.tsv
done
# Make list of batches
for i in $( seq 1 ${n_splits} ); do
  echo -e "gnomAD-SV_v3.wgd_ploidy_batch.batch_${i}"
done > ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list


### Prepare inputs for module 00b for all batches
# Make zip with dependencies
if [ -e ${WRKDIR}/dependencies.zip ]; then
  rm ${WRKDIR}/dependencies.zip
fi
cd ${GATKSV}/module00b/ \
&& zip dependencies.zip *wdl \
&& mv dependencies.zip ${WRKDIR}/ \
&& cd -
# Make template json for inputs
cd ${GATKSV}/module00b/ \
&& womtool inputs Module00b.wdl \
> ${WRKDIR}/json/template.input.json \
&& cd -
# Write input json for each batch
while read batch; do
  manifest=${WRKDIR}/batches/${batch}.samples_wCovPath.tsv
  samples=$( awk -v ORS="," '{ print "\""$1"\"" }' $manifest \
             | sed -e 's/^/\[/g' -e 's/,$/\]/g' )
  covfiles=$( awk -v ORS="," '{ print "\""$2"\"" }' $manifest \
              | sed -e 's/^/\[/g' -e 's/,$/\]/g' )
  cat <<EOF > ${WRKDIR}/json/${batch}.input.json
{
  "Module00b.samples": ${samples},
  "Module00b.counts": ${covfiles},
  "Module00b.run_vcf_qc": false,
  "Module00b.sv_pipeline_docker": "rlcollins/sv-pipeline:wgd_sample_id_hotfix-d12bcf",
  "Module00b.genome_file": "gs://talkowski-sv-gnomad-resources/lists/hg38.genome",
  "Module00b.sv_mini_docker": "rlcollins/sv-base-mini:wgd_sample_id_hotfix-d12bcf",
  "Module00b.wgd_scoring_mask": "gs://talkowski-sv-gnomad-resources/rlc_refs/wgd_scoring_mask.hg38.gnomad_v3.bed",
  "Module00b.run_ploidy": true,
  "Module00b.sv_pipeline_qc_docker": "rlcollins/sv-pipeline-qc:wgd_sample_id_hotfix-d12bcf",
  "Module00b.batch": "${batch}"
}
EOF
done < ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list


### Run first 20 batches (as a test)
# Submit
while read batch; do
  cromshell -t 20 submit \
    ${GATKSV}/module00b/Module00b.wdl \
    ${WRKDIR}/json/${batch}.input.json \
    ${cromwell_options} \
    ${WRKDIR}/dependencies.zip \
  > ${WRKDIR}/logs/${batch}.cromshell_submission.log
done < <( head -n20 ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list )
# Check status
while read batch; do
  cromshell -t 20 metadata \
    $( cat ${WRKDIR}/logs/${batch}.cromshell_submission.log \
       | sed '1d' | jq .id | tr -d '"' ) \
  | jq .status
done < <( head -n20 ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list )
# Download ploidy & WGD results for first test batch (to confirm workflow outputs)
batch=$( head -n1 ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list )
cromshell -t 20 metadata \
  $( cat ${WRKDIR}/logs/${batch}.cromshell_submission.log \
     | sed '1d' | jq .id | tr -d '"' ) \
| jq .outputs | grep 'ploidy_plots\|WGD_dist' \
| awk '{ print $2 }' | tr -d '"' | sed 's/,$//g' \
| gsutil -m cp -I ~/scratch/
# Download and process first 5 batches as confirmation
head -n5 ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list \
> ~/scratch/test_batches.list
$BASEDIR/gnomad_sv_v3/sample_qc_and_batching/module00b_wgd_ploidy_rerun_jan2020/postprocess_wgd_ploidy.sh \
  ~/scratch/test_batches.list \
  ~/scratch/test_wgd_ploidy_metadata.merged.tsv.gz


### Run all batches
# Submit
while read batch; do
  cromshell -t 20 submit \
    ${GATKSV}/module00b/Module00b.wdl \
    ${WRKDIR}/json/${batch}.input.json \
    ${cromwell_options} \
    ${WRKDIR}/dependencies.zip \
  > ${WRKDIR}/logs/${batch}.cromshell_submission.log
done < ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list
# Check status
while read batch; do
  cromshell -t 20 metadata \
    $( cat ${WRKDIR}/logs/${batch}.cromshell_submission.log \
       | sed '1d' | jq .id | tr -d '"' ) \
  | jq .status
done < ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list \
> ~/scratch/cromshell_status_check.txt
# Download and process all batches
$BASEDIR/gnomad_sv_v3/sample_qc_and_batching/module00b_wgd_ploidy_rerun_jan2020/postprocess_wgd_ploidy.sh \
  ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list \
  ${BASEDIR}/data/mod00b_wgd_ploidy/gnomAD-SV_v3.wgd_ploidy_metadata.tsv.gz \
> ~/scratch/postprocess_wgd_ploidy.stdout.log


### May 2020: Rerun of 606 samples with truncated IDs due to periods in sample IDs
xz_outbucket="gsutil -m ls gs://talkowski-sv-gnomad-output/zero/Module00b/8ffe3739-dfa8-4b3c-a947-37048f0bf16b/"
gs_ploidy_matrix="${xz_outbucket}**gnomAD-SV_v3_sampleid_adjustment_ploidy_matrix.bed.gz"
gs_wgd_scores="${xz_outbucket}**gnomAD-SV_v3_sampleid_adjustment_WGD_scores.txt.gz"
if ! [ -e ${WRKDIR}/makeup606_april2020 ]; then
  mkdir ${WRKDIR}/makeup606_april2020
fi
gsutil -m cp \
  ${gs_ploidy_matrix} \
  ${gs_wgd_scores} \
  ${WRKDIR}/makeup606_april2020/


### Also May 2020: Addition of 698 1kGB relatives
if ! [ -e ${WRKDIR}/1kGP_698_relatives_may2020 ]; then
  mkdir ${WRKDIR}/1kGP_698_relatives_may2020
fi
if [ -e ${WRKDIR}/1kGP_698_relatives_may2020/output_shards ]; then
  rm -rf ${WRKDIR}/1kGP_698_relatives_may2020/output_shards
fi
mkdir ${WRKDIR}/1kGP_698_relatives_may2020/output_shards
# Download WGD scores
gsutil -m cp \
  gs://talkowski-sv-gnomad-output/1KGP/Module00b/WGD/1KGP_698.**_WGD_scores.txt.gz \
  ${WRKDIR}/1kGP_698_relatives_may2020/output_shards/
# Download ploidy results
gsutil -m cp \
  gs://talkowski-sv-gnomad-output/1KGP/Module00b/Ploidy/1KGP_698.**_ploidy_plots.tar.gz \
  gs://talkowski-sv-gnomad-output/1KGP/Module00b/Ploidy/1KGP_698.**_ploidy_matrix.bed.gz \
  ${WRKDIR}/1kGP_698_relatives_may2020/output_shards/
# Postprocess all data
$BASEDIR/gnomad_sv_v3/sample_qc_and_batching/module00b_wgd_ploidy_rerun_jan2020/postprocess_wgd_ploidy.1kGB_698_relatives.may2020.sh