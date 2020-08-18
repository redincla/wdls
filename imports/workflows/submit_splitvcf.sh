#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=SplitVCF
#SBATCH -t 00-01:00
#SBATCH --mem-per-cpu=10G
#SBATCH -o slurm.%N.%j.log

#sample_ID=$1
module add Development/java/1.8.0_232
module add R/3.6.1

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
#     -V /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/sample/old_ACorig/0620.filtered.annotated_CardioGenes_${sample_ID}.vcf.gz \
#     -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
#     -F 1000g2015aug_all -F AAChange.ensGene -F AAChange.refGene -F AC_Orig -F AF_Orig -F AF_afr -F AF_ami -F AF_amr -F AF_asj -F AF_eas -F AF_female -F AF_fin -F AF_male -F AF_nfe -F AF_popmax -F AN_Orig -F AS_FS -F AS_MQ \
#     -F AS_QD -F CADD_phred -F CLINSIG -F Confidence -F DP -F ExonicFunc.ensGene -F ExonicFunc.refGene -F Func.ensGene -F Func.refGene -F GERP++_RS -F Gene.ensGene -F Gene.refGene -F Inheritance \
#     -F Kaviar_AF -F MutationAssessor_pred -F MutationAssessor_score -F MutationTaster_pred -F MutationTaster_score -F NEGATIVE_TRAIN_SITE -F POSITIVE_TRAIN_SITE -F Polyphen2_HDIV_pred -F Polyphen2_HDIV_score \
#     -F Polyphen2_HVAR_pred -F Polyphen2_HVAR -F QD -F SIFT_pred -F SIFT_score -F n.PLP.ClinVar -F n.PLP.LoF.ClinVar -F n.PLP.mis.ClinVar -F non_cancer_AF_popmax -F oe.LoF.upper \
#     -F pLI -F phastCons100way_vertebrate -F phyloP100way_vertebrate -F phyloP20way_mammalian \
#     -O /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/sample/old_ACorig/0620.filtered.annotated_CardioGenes_${sample_ID}.vcf.tsv

java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
VariantsToTable \
     -V /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Delly/0620.delly_CardioCRAM.vcf.gz \
     -F CHROM -F POS -F ID -F REF ALT -F QUAL -F FILTER -F CIEND -F CIPOS -F CHR2 -F END -F PE \
     -F MAPQ -F SR -F SRQ -F CONSENSUS -F CT -F IMPRECISE -F PRECISE -F SVTYPE -F SVMETHOD -F INSLEN -F HOMLEN \
     -GF GT -GF GL -GF GQ -GF FT -GF RC -GF RCL -GF RCR \
     -GF CN -GF DR -GF DV -GF RR -GF RV \
     -O /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Delly/0620.delly_CardioCRAM.vcf.tsv

# awk 'NR==1; NR > 1{s=0; for (i=2;i<=NF;i++) s+=$i; if ($11 < 4) print }' 0620.filtered.annotated_CardioGenes_H3_Plate_103.vcf.tsv > 0620.filtered.annotated_CardioGenes_H3_Plate_103_AC-orig-filtered.vcf.tsv


#java -Xmx8g -jar /software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar \
#SelectVariants \
#-R /data/chuv/LABO/redin/references/hg38/Homo_sapiens_assembly38.fasta \
#-V /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/sample/0620.filtered.annotated_CardioGenes_${sample_ID}.tmp.vcf.gz \
#-O /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/sample/0620.filtered.annotated_CardioGenes_${sample_ID}.vcf.gz \
#-select "AC_Orig < 5" 

#rm /home/credin/scratch/WGS/data_and_refs/data/Raw_FQ/0620/Annotation/sample/0620.filtered.annotated_CardioGenes_${sample_ID}.tmp.vcf.gz