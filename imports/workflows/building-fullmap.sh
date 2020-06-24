#!/bin/bash

########################
#    Creating full.map  #
########################
# col 1: fastq_1 , col 2: fastq_2 , col3: RG, col4: lib ID, col 5: PU, col6: run date, col7: platform, col8: seq center

export BASEDIR=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ
export WRKDIR=${BASEDIR}/0620

# checks for appropriate input
if [ $# -eq 2 ]; then
 sample_list=$1 #col1: sample list, col2: lib preparation date
 FQDIR=$2 #full directory of FQ files
else
 exit 1
fi

cd ${FQDIR}
num_col=$( echo "${FQDIR}" | awk -F'/' '{print NF; exit}')
run_date=$( echo "${FQDIR}" | cut -f ${num_col} -d/ | cut -f 1 -d_)
FQDIR=$( echo "${FQDIR}" | sed "s+/$++g")  #replace last '/' if present in given fqdir path in order not to mess up the final fq absolute paths

while read sample_ID lib_date; do  ## build sample map
ls -l | awk '{print $9}' > tmp
grep -P "${sample_ID}.*R1" tmp > ${sample_ID}.fq1_list #get list of fq_R1
   while read fq_R1; do
        fq_R2=$( echo "${fq_R1}" | sed "s+R1+R2+g" )
        RG=$( echo "${fq_R1}" | sed "s+_R1+/+g" | cut -f 1 -d/ )
        PU=$( zcat ${fq_R1} | cut -f 3,4,10 -d: | head -n 1 )
        printf "%s/%s\t%s/%s\t%s\t%s\t%s\t%s\tillumina\tH2030GC\n" "${FQDIR}" "${fq_R1}" "${FQDIR}" "${fq_R2}" "${RG}" "${lib_date}" "${PU}" "${run_date}" >> ${FQDIR}/${sample_ID}.full_map.tsv
    done < ${sample_ID}.fq1_list
    rm ${sample_ID}.fq1_list
    rm tmp
    
    if [ -f "${BASEDIR}/${sample_ID}.full_map.tsv" ]; then
        cat ${BASEDIR}/${sample_ID}.full_map.tsv ${FQDIR}/${sample_ID}.full_map.tsv > tmp2
        mv tmp2 ${BASEDIR}/${sample_ID}.full_map.tsv
        rm ${FQDIR}/${sample_ID}.full_map.tsv
    else
        mv ${FQDIR}/${sample_ID}.full_map.tsv ${BASEDIR}
    fi
done < ${sample_list}