while read gene_symbol gene_mim disease_name disease_mim DDD_category allelic_requirement mutation_consequence phenotypes organ_specificity_list pmids panel prev_symbols hgnc_id gene_disease_pair_entry_date
do
    while read chr start end Gene
    do
        if grep -qi "${GeneSymbol}" DDG2P_clean_sorted.bed
        then
            if [ "${gene_symbol}" == "${Gene}" ]
            then
                printf "${chr}\t${start}\t${end}\t${Gene}\t${gene_symbol}\t${gene_mim}\t${disease_name}\t${disease_mim}\t${DDD_category}\t${allelic_requirement}\t${mutation_consequence}\t${organ_specificity_list}\t${gene_disease_pair_entry_date}\n" >> DDG2P_010920.FULL.bed
            fi
        else
        printf "${GeneSymbol}\n"
        fi
    done < DDG2P_clean_sorted.bed
done < DDG2P_1_9_2020_sorted.tsv


vcfanno ~/scratch/WGS/wdls/imports/config-files/vcfanno-config.toml  /scratch/beegfs/PRTNR/CHUV/chuv_medp/WGS/wdls/cromwell-executions/JointAnnotation/7494e2ff-c557-45ec-bddf-c4042e183c7d/call-vcfannoScatteredVCF/shard-0/inputs/-1725729894/0620.filtered.0.hg38_multianno.vcf.gz > "test2.hg38_vcfanno.vcf"