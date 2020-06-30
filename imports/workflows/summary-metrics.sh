#!/bin/bash

########################
#   gather WGS metrics   #
########################

# Copyright (c) 2020 Claire Redin and The CHUV Precision Medicine Unit, Lausanne, Switzerland.
# Contact: Claire Redin <claire.redin@chuv.ch>

# checks for appropriate input
if [ $# -eq 1 ]; then
 sample_list=$1 #full path to file: list of samples to gather WGS metrics on

else
 echo -e "\n\nGather WGS summary metrics by sample\n\nAuthor: Claire Redin (claire.redin@chuv.ch)\n\n"
 echo "Usage:"
 echo "summary-metrics.sh [sample_list] "
 echo "sample_list: full path to file: list of samples to gather WGS metrics on"
 exit 1
fi

export BASEDIR=/home/credin/scratch/WGS/data_and_refs/data/Raw_FQ
export WRKDIR=${BASEDIR}/0620

while read sample_ID; do
    cd ${WRKDIR}
################ CROMWELL #############
## Check that workflow succeed
    if grep -q "WorkflowSucceededState" slurm.${sample_ID}.*.log ; then
        printf "%s: SUCCESS\n" "${sample_ID}">> ${WRKDIR}/Agregate_workflow.log
        CONTAMINATION=$(grep "WholeGenomeGermlineSingleSample.contamination" slurm.${sample_ID}.*.log | cut -f 2 -d: | cut -f 1 -d, | head -n 1)
        cd ${WRKDIR}/${sample_ID}/Processed

################ PER UBAM #############
## Check that all ubam file have no errors
#    ls -l | awk '{print $9}' | grep .bam.validation_report > ubam.validation.list
#    while read report; do
#        if ! grep -q "No errors found" ${report}; then
#            echo -e "${sample_ID}: errors in ${report}\n" >> ${WRKDIR}/Agregate_bam.log
#        fi
#    done < ${ubam.validation.list}
    #rm ubam.validation.list

################ PER SAMPLE #############
## gather summary metrics per sample

        if ! grep -q "No errors found" *cram.validation_report; then
            echo -e "${sample_ID}: errors in cram file\n" >> ${WRKDIR}/Agregate_bam.log
        fi

###duplicate metrics
        header_line=$(grep -n "UNPAIRED_READS_EXAMINED" ${sample_ID}.duplicate_metrics | cut -f1 -d: )
        num_line=$(( ${header_line} + 1))
        READ_PAIRS_EXAMINED=$(sed -n ${num_line}p ${sample_ID}.duplicate_metrics | cut -f 3)
        PERCENT_DUPLICATION=$(sed -n ${num_line}p ${sample_ID}.duplicate_metrics | cut -f 9)

###alignment metrics
        num_line=$(grep -Pn "^PAIR\t" ${sample_ID}.alignment_summary_metrics | cut -f1 -d: | head -1 )
        TOTAL_READS=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 2)
        PF_READS=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 3)
        PCT_PF_READS=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 4)
        PCT_PF_READS_ALIGNED=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 7)
        PF_HQ_ALIGNED_READS=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 9)
        PCT_PF_HQ_ALIGNED_READS=$(( ${PF_HQ_ALIGNED_READS} / ${PF_READS} ))  
        PF_HQ_ERROR_RATE=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 14)
        PF_INDEL_RATE=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 15)
        MEAN_READ_LENGTH=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 16)
        PCT_READS_ALIGNED_IN_PAIRS=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 18)
        PCT_PF_READS_IMPROPER_PAIRS=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 20)
        BAD_CYCLES=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 21)
        STRAND_BALANCE=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 22)
        PCT_CHIMERAS=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 23)
        PCT_ADAPTER=$(sed -n ${num_line}p ${sample_ID}.alignment_summary_metrics | cut -f 24)

###insert size metrics
        header_line=$(grep -n "MEDIAN_INSERT_SIZE" ${sample_ID}.insert_size_metrics | cut -f1 -d: )
        num_line=$(( ${header_line} + 1))
        MEDIAN_INSERT_SIZE=$(sed -n ${num_line}p ${sample_ID}.insert_size_metrics | cut -f 1)
        MEDIAN_ABSOLUTE_DEVIATION=$(sed -n ${num_line}p ${sample_ID}.insert_size_metrics | cut -f 3)

###wgs metrics
        header_line=$(grep -n "GENOME_TERRITORY" ${sample_ID}.wgs_metrics | cut -f1 -d: )
        num_line=$(( ${header_line} + 1))
        MEDIAN_COVERAGE=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 4)
        MAD_COVERAGE=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 5)
        PCT_EXC_ADAPTER=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 6)
        PCT_EXC_MAPQ=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 7)
        PCT_EXC_DUPE=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 8)
        PCT_EXC_UNPAIRED=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 9)
        PCT_EXC_BASEQ=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 10)
        PCT_EXC_TOTAL=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 13)
        PCT_10X=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 16)
        PCT_20X=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 18)
        PCT_30X=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 20)
        HET_SNP_SENSITIVITY=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 28)
        HET_SNP_Q=$(sed -n ${num_line}p ${sample_ID}.wgs_metrics | cut -f 29)

###variant calling metrics
        header_line=$(grep -n "SAMPLE_ALIAS" ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f1 -d: )
        num_line=$(( ${header_line} + 1))
        HET_HOMVAR_RATIO=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 2)
        PCT_GQ0_VARIANTS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 3)
        TOTAL_SNPS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 6)
        NUM_IN_DB_SNP=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 7)
        NOVEL_SNPS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 8)
        PCT_DBSNP=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 10)
        DBSNP_TITV=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 11)
        NOVEL_TITV=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 12)
        TOTAL_INDELS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 13)
        NOVEL_INDELS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 14)
        PCT_DBSNP_INDELS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 16)
        DBSNP_INS_DEL_RATIO=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 18)
        NOVEL_INS_DEL_RATIO=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 19)
        TOTAL_MULTIALLELIC_SNPS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 20)
        NUM_IN_DB_SNP_MULTIALLELIC=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 21)
        TOTAL_COMPLEX_INDELS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 22)
        NUM_IN_DB_SNP_COMPLEX_INDELS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 23)
        NUM_SINGLETONS=$(sed -n ${num_line}p ${sample_ID}_PLUMBING.variant_calling_detail_metrics | cut -f 25)

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"  \
        "${sample_ID}" "${CONTAMINATION}" "${TOTAL_READS}" "${PF_READS}" "${PCT_PF_READS}" "${PCT_PF_READS_ALIGNED}" "${PCT_PF_HQ_ALIGNED_READS}" "${PF_HQ_ERROR_RATE}" "${PF_INDEL_RATE}" "${MEAN_READ_LENGTH}" "${PCT_READS_ALIGNED_IN_PAIRS}" "${PCT_PF_READS_IMPROPER_PAIRS}" "${BAD_CYCLES}" "${STRAND_BALANCE}" "${PCT_CHIMERAS}" "${PCT_ADAPTER}" \
        "${PERCENT_DUPLICATION}" \
        "${MEDIAN_INSERT_SIZE}" "${MEDIAN_ABSOLUTE_DEVIATION}" \
        "${MEDIAN_COVERAGE}" "${MAD_COVERAGE}" "${PCT_EXC_ADAPTER}" "${PCT_EXC_MAPQ}" "${PCT_EXC_DUPE}" "${PCT_EXC_UNPAIRED}" "${PCT_EXC_BASEQ}" "${PCT_EXC_TOTAL}" "${PCT_10X}" "${PCT_20X}" "${PCT_30X}" "${HET_SNP_SENSITIVITY}" "${HET_SNP_Q}" \
        "${HET_HOMVAR_RATIO}" "${PCT_GQ0_VARIANTS}" "${TOTAL_SNPS}" "${NOVEL_SNPS}" "${PCT_DBSNP}" "${DBSNP_TITV}" "${NOVEL_TITV}" "${TOTAL_INDELS}" "${NOVEL_INDELS}" "${PCT_DBSNP_INDELS}" "${DBSNP_INS_DEL_RATIO}" "${NOVEL_INS_DEL_RATIO}" "${TOTAL_MULTIALLELIC_SNPS}" "${NUM_IN_DB_SNP_MULTIALLELIC}" "${TOTAL_COMPLEX_INDELS}" "${NUM_IN_DB_SNP_COMPLEX_INDELS}" "${NUM_SINGLETONS}" >> ${WRKDIR}/Agregate_WGS_metrics.tsv

    else
        printf "%s: Workflow Failed\n" "${sample_ID}" >> ${WRKDIR}/Agregate_workflow.log
    fi

done < ${sample_list}