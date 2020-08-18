#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=rliftover
#SBATCH -t 1-00:00
#SBATCH --mem-per-cpu=2G

#module add Utility/UCSC-utils/359

#grep "^#" hg19_spidex.txt > header
#grep -v '#' hg19_spidex.txt | awk -v s=1 '{print "chr"$1"\t"$2"\t"$3+s"\t"$4"\t"$5"\t"$6"\t"$7;}' > hg19_spidex_tmp.txt

#liftOver /home/credin/refs/tools/annovar/AnnovarDB/humandb/hg38/hg19_spidex_tmp.txt \
#~/refs/references/hg38/hg19ToHg38.over.chain.gz \
#/home/credin/refs/tools/annovar/AnnovarDB/humandb/hg38/hg19_spidex_tmp_liftedhg38.txt \
#-bedPlus=3 unMapped

#cat $header | awk -v s=1 '{print $1"\t"$2"\t"$3-s"\t"$4"\t"$5"\t"$6"\t"$7;}' hg19_spidex_tmp_liftedhg38.txt | sed 's@^chr@@g' > hg38_spidex_lifted-from-hg19.txt

sed 's@^chr@@g' /home/credin/refs/tools/annovar/AnnovarDB/humandb/hg38/hg38_gnomad30_genome.txt > /home/credin/refs/tools/annovar/AnnovarDB/humandb/hg38/hg38_gnomad30_genome_tmp.txt

#rm header unMapped hg19_spidex_tmp.txt hg19_spidex_tmp_liftedhg38.txt