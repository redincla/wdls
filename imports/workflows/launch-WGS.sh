#!/usr/bin/env bash

########################
#    WGS pipeline v1    #
########################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit
# Contact: Claire Redin <claire.redin@chuv.ch>

### Set local parameters
export BASEDIR=/scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/data_and_refs/data/Raw_FQ/
export WRKDIR=${BASEDIR}/${cohort_name}_0620

if ! [ -e $WRKDIR ]; then
	mkdir $WRKDIR
fi

cd ${WRKDIR}

#export cromwell_options="${BASEDIR}/cromwell_options.json"

### Prepare folder structure
while read sample_ID; do
    if ! [ -e $sample_ID ]; then
	    mkdir $sample_ID
     if ! [ -e ${sample_ID}/QC ]; then
	    mkdir ${sample_ID}/QC
     if ! [ -e ${sample_ID}/Raw ]; then
	    mkdir ${sample_ID}/Raw
     if ! [ -e ${sample_ID}/Processed ]; then
	    mkdir ${sample_ID}/Processed

    mv ${BASEDIR}/${sample_ID}*fastq.gz ${WRKDIR}/${sample_ID}/Raw/
done < (cut -f 1 ${full_map} | sort -u )

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
  "ConvertPairedFastQsToUnmappedBamWf.GATK": "/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar",
  "ConvertPairedFastQsToUnmappedBamWf.PICARD": "/software/UHTS/Analysis/picard-tools/2.21.8/bin/picard.jar",
  "ConvertPairedFastQsToUnmappedBamWf.full_map": ${full_map},
  "ConvertPairedFastQsToUnmappedBamWf.make_fofn": true,
  "ConvertPairedFastQsToUnmappedBamWf.cohort_name": ${cohort_name}
}
EOF
done < ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list


### Run first 20 batches (as a test)
# Submit

launch Convert scripts
launch WGS script 

while read batch; do
  cromshell -t 20 submit \
    ${GATKSV}/module00b/Module00b.wdl \
    ${WRKDIR}/json/${batch}.input.json \
    ${cromwell_options} \
    ${WRKDIR}/dependencies.zip \
  > ${WRKDIR}/logs/${batch}.cromshell_submission.log
done < <( head -n20 ${WRKDIR}/gnomAD-SV_v3.wgd_ploidy_batches.list )


    place QC output here
    place all final (not intermediate) processed-files here
