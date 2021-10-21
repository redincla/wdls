#!/bin/bash

########################
#    Creating full.map  #
########################
# col 1: fastq_1 , col 2: fastq_2 , col3: RG, col4: lib ID, col 5: PU, col6: run date, col7: platform, col8: seq center
# rename fq files including seq date? So that when copying to Raw directory later on there are not overwriting risks...

# checks for appropriate input
if [ $# -eq 3 ]; then
# sample_list=$1 #col1: sample list, deprecated: col2: lib preparation date
 BASEDIR=$1 #full directory of FQ files
 DESTDIR=$2 #full directory to drop sample_maps
 run_dates=$3 #file with correspondance between: col1: run_ID, col2: run_date
else
 echo -e "\n\nBuilding sample map script by sample\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "building-fullmap.sh [fq directory] [output directory] [run_dates]"
 echo "- fq directory: full path"
 echo "- output directory: full path"
 echo "- run_dates: 2-columns file, giving correspondance between: col1: run_ID, col2: run_date"
 exit 1
fi

BASEDIR=$( echo "${BASEDIR}" | sed "s+/$++g") #replace last '/' if present in given fqdir path in order not to mess up the final fq absolute paths

export DESTDIR
export BASEDIR

cd ${BASEDIR}

####gathering all necessary files to reorganize fq files by run_ID folders
#get list of fqs
find . -name "*.gz" | cut -f 2 -d/ > fq_list

# get list of run_IDs and corresponding run_ID for each fq for later reorganization
while read fq; do
run_ID=$(zcat ${fq} | cut -f 1,2 -d: | head -n 1 | sed "s+:+.+g" | sed "s+@++g") 
zcat ${fq} | cut -f 1,2 -d: | head -n 1 | sed "s+:+.+g" | sed "s+@++g" >> run_ID_list
printf "%s\t%s\n" "${fq}" "${run_ID}" >> fq_dest-folder
done < fq_list

sort -u run_ID_list > run_ID_list_uniq
mv run_ID_list_uniq run_ID_list

#create list of uniq run_ID folders
while read line; do
    if ! [ -e ${line} ]; then
	    mkdir ${line}
    fi
done < run_ID_list

#move fq into corresponding run_ID destination folder
while read fq run_ID; do
    mv ${fq} ${run_ID}
done < fq_dest-folder


while read runID; do

    FQDIR="${BASEDIR}/${runID}"
    cd ${FQDIR}
    FQDIR=$( echo "${FQDIR}" | sed "s+/$++g")  #replace last '/' if present in given fqdir path in order not to mess up the final fq absolute paths
    
    #get list of sample ID
    find . -name "*.gz"  | awk -F_NGS '{print $1 }' | cut -f 2 -d/ | sort -u > sample_list

    while read sample_ID; do  ## build sample map
        find . -name "${sample_ID}_*_R1_*.fastq.gz" | sed "s+^./++" > ${sample_ID}.fq1_list #get list of fq_R1
        while read fq_R1; do
            fq_R2=$( echo "${fq_R1}" | sed "s+R1+R2+g" )
            run_ID=$( zcat ${fq_R1} | cut -f 1,2 -d: | head -n 1 | sed "s+:+.+g" | sed "s+@++g" ) #/!\ has to be NNNNNN date format
            lib_index=$(( $(echo "${fq_R1}" | awk '{print index($0,"NGS")}') - 1 ))
            lib_ID=$( echo "${fq_R1:lib_index}" |  cut -f 1 -d_ )
            RG=$( echo ${sample_ID}.$(zcat ${fq_R1} | cut -f 3,4 -d: | head -n 1 | sed "s+:+.+g") )
            PU=$( zcat ${fq_R1} | cut -f 3,4,10 -d: | head -n 1 | sed "s+:+.+g" )
            printf "%s/%s_%s\t%s/%s_%s\t%s\t%s\t%s\t%s\tillumina\tH2030GC\n" "${FQDIR}" "${run_ID}" "${fq_R1}" "${FQDIR}" "${run_ID}" "${fq_R2}" "${RG}" "${lib_ID}" "${PU}" "${run_ID}" >> ${FQDIR}/${sample_ID}.full_map.tsv
#           printf "%s/%s\t%s/%s\t%s\t%s\t%s\t%s\tillumina\tH2030GC\n" "${FQDIR}" "${fq_R1}" "${FQDIR}" "${fq_R2}" "${RG}" "${lib_ID}" "${PU}" "${run_ID}" >> ${FQDIR}/${sample_ID}.full_map.tsv
            mv ${fq_R1} ${run_ID}_${fq_R1} #to avoid FQ name duplicates when coming from different run IDs
            mv ${fq_R2} ${run_ID}_${fq_R2}
        done < ${sample_ID}.fq1_list
        rm ${sample_ID}.fq1_list

        while read run date; do
            sed -i "s+\t${run}\t+\t${date}\t+g" ${FQDIR}/${sample_ID}.full_map.tsv
        done < ${run_dates}
    
        if [ -f "${DESTDIR}/${sample_ID}.full_map.tsv" ]; then
            cat ${DESTDIR}/${sample_ID}.full_map.tsv ${FQDIR}/${sample_ID}.full_map.tsv > tmp2
            mv tmp2 ${DESTDIR}/${sample_ID}.full_map.tsv
            rm ${FQDIR}/${sample_ID}.full_map.tsv
        else
            mv ${FQDIR}/${sample_ID}.full_map.tsv ${DESTDIR}
        fi
    done < sample_list

done < run_ID_list