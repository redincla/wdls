#!/bin/bash

########################
#    Creating full.map  #
########################
# col 1: fastq_1 , col 2: fastq_2 , col3: RG, col4: lib ID, col 5: PU, col6: run date, col7: platform, col8: seq center
# rename fq files including seq date? So that when copying to Raw directory later on there are not overwriting risks...

#export BASEDIR=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/TCAG
#export WRKDIR=${BASEDIR}/0920

# checks for appropriate input
if [ $# -eq 3 ]; then
 sample_list=$1 #col1: sample list, deprecated: col2: lib preparation date
 FQDIR=$2 #full directory of FQ files
 BASEDIR=$3 #full directory to drop sample_maps
else
 echo -e "\n\nBuilding sample map script by sample\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "building-fullmap.sh [sample_list] [fq directory] [output directory]"
 echo "sample_list: Column 1= list of samples to process"
 echo "fq directory: full path"
 echo "output directory: full path"
 exit 1
fi

export BASEDIR
cd ${FQDIR}
#num_col=$( echo "${FQDIR}" | awk -F'/' '{print NF; exit}')
#run_date=$( echo "${FQDIR}" | cut -f ${num_col} -d/ | cut -f 1 -d_)
FQDIR=$( echo "${FQDIR}" | sed "s+/$++g")  #replace last '/' if present in given fqdir path in order not to mess up the final fq absolute paths

while read sample_ID; do  ## build sample map
ls -l | awk '{print $9}' > tmp
grep -P "${sample_ID}_.*R1" tmp > ${sample_ID}.fq1_list #get list of fq_R1
   while read fq_R1; do
        fq_R2=$( echo "${fq_R1}" | sed "s+R1+R2+g" )
        run_ID=$( zcat ${fq_R1} | cut -f 1,2 -d: | head -n 1 | sed "s+:+.+g" | sed "s+@++g" ) #/!\ has to be NNNNNN date format
        lib_ID=$( echo "${fq_R1}" |  cut -f 2 -d_ )
        RG=$( zcat ${fq_R1} | cut -f 3,4 -d: | head -n 1 | sed "s+:+.+g" )
        PU=$( zcat ${fq_R1} | cut -f 3,4,10 -d: | head -n 1 | sed "s+:+.+g" )
        echo "${FQDIR} ${fq_R1} ${FQDIR} ${fq_R2} ${RG} ${lib_ID} ${PU} ${run_ID}"
        printf "%s/%s_%s\t%s/%s_%s\t%s\t%s\t%s\t%s\tillumina\tH2030GC\n" "${FQDIR}" "${run_ID}" "${fq_R1}" "${FQDIR}" "${run_ID}" "${fq_R2}" "${RG}" "${lib_ID}" "${PU}" "${run_ID}" >> ${FQDIR}/${sample_ID}.full_map.tsv
#        printf "%s/%s\t%s/%s\t%s\t%s\t%s\t%s\tillumina\tH2030GC\n" "${FQDIR}" "${fq_R1}" "${FQDIR}" "${fq_R2}" "${RG}" "${lib_ID}" "${PU}" "${run_ID}" >> ${FQDIR}/${sample_ID}.full_map.tsv
        mv ${fq_R1} ${run_ID}_${fq_R1}
        mv ${fq_R2} ${run_ID}_${fq_R2}
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