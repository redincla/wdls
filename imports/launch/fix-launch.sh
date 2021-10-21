#!/bin/bash

########################
#    WGS pipeline v1    #
########################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit, Lausanne, Switzerland.
# Contact: Claire Redin <claire.redin@chuv.ch>

# checks for appropriate input
if [ $# -eq 1 ]; then
 sample_list=$1 #list of samples to reprocess, one sample / line
else
 echo "Usage:\n"
 echo "fix-launch.sh [sample_list] "
 echo "sample_list: full path to list of samples to process "
 exit 1
fi

### Set local parameters
export BASEDIR=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ
export WRKDIR=${BASEDIR}/0620

cd ${WRKDIR}

while read sample_ID; do
    mv ${WRKDIR}/${sample_ID}/Raw/${sample_ID}*fastq.gz  ${BASEDIR}/
#    mv ${WRKDIR}/${sample_ID}/Raw/${sample_ID}.full_map.tsv ${BASEDIR}/
#    sed -i "s+${BASEDIR}//0620/${sample_ID}/Raw+${BASEDIR}/+g" ${BASEDIR}/${sample_ID}.full_map.tsv
    cd ${WRKDIR}
    rm -r $sample_ID
done < ${sample_list}