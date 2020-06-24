#!/bin/bash

########################
#    WGS pipeline v1    #
########################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit, Lausanne, Switzerland.
# Contact: Claire Redin <claire.redin@chuv.ch>

# checks for appropriate input
if [ $# -eq 1 ]; then
 sample_list=$1 #full path to list of samples to process from same batch, one sample / line
# cohort_name=$2 #batch name - should be the same as folder name under which FQ.gz files are stored
else
 echo -e "\n\nLaunching WGS parallelizer script by sample\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "launch-WGS.sh [sample_list] "
 echo "sample_list: full path to list of samples to process "
 exit 1
fi

### Set local parameters
export BASEDIR=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ
export WRKDIR=${BASEDIR}/0620

#cd ${WRKDIR}

#export cromwell_options="${BASEDIR}/cromwell_options.json"

### Prepare folder structure  -> Adapt full map with relative paths, not absolute paths since FQ files will be relocated
while read sample_ID; do
#    mv ${WRKDIR}/${sample_ID}/Raw/${sample_ID}*fastq.gz  ${BASEDIR}/
#    mv ${WRKDIR}/${sample_ID}/Raw/${sample_ID}.full_map.tsv ${BASEDIR}/
#    sed -i "s+${BASEDIR}//0620/${sample_ID}/Raw+${BASEDIR}/+g" ${BASEDIR}/${sample_ID}.full_map.tsv
    cd ${WRKDIR}
    rm -r $sample_ID
done < ${sample_list}