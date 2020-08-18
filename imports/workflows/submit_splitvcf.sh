#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=SplitVCF
#SBATCH -t 00-04:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o slurm.%N.%j.log

#sample_ID=$1
module add Development/java/1.8.0_232

java -Dconfig.file=/home/credin/.cromwell.conf_new -jar /software/Utility/cromwell/47/bin/cromwell-47.jar run /home/credin/scratch/WGS/wdls/imports/workflows/FilterVCF.wdl -i /home/credin/scratch/WGS/wdls/imports/config-files/FilterVCF.inputs.json -o /home/credin/scratch/WGS/wdls/imports/config-files/FilterVCF.options.json

#java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
#SelectVariants \
#--intervals /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/list.bed \
#-R /data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta \
#-V /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/0620.filtered.annotated.vcf.gz \
#-O /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/sample/0620.filtered.annotated_CardioGenes_${sample_ID}.tmp.vcf.gz \
#--exclude-non-variants \
#--keep-original-ac \
#-sn ${sample_ID} 

#java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
#VariantsToTable \
#    -V /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Delly/0620.delly_CardioCRAM.vcf.gz \
#     -F CHROM -F POS -F ID -F REF ALT -F QUAL -F FILTER -F CIEND -F CIPOS -F CHR2 -F END -F PE \
#     -F MAPQ -F SR -F SRQ -F CONSENSUS -F CT -F IMPRECISE -F PRECISE -F SVTYPE -F SVMETHOD -F INSLEN -F HOMLEN \
#     -GF GT -GF GL -GF GQ -GF FT -GF RC -GF RCL -GF RCR \
#     -GF CN -GF DR -GF DV -GF RR -GF RV \
#     -O /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Delly/0620.delly_CardioCRAM.vcf.tsv

# awk 'NR==1; NR > 1{s=0; for (i=2;i<=NF;i++) s+=$i; if ($11 < 4) print }' 0620.filtered.annotated_CardioGenes_H3_Plate_103.vcf.tsv > 0620.filtered.annotated_CardioGenes_H3_Plate_103_AC-orig-filtered.vcf.tsv

#rm /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/sample/0620.filtered.annotated_CardioGenes_${sample_ID}.tmp.vcf.gz

###step 4: filter on cardiogenes