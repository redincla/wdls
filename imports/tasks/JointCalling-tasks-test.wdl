version 1.0

############
### check that there are no sample doublons
############
task CheckSamplesUnique {
  input {
    File sample_name_map
    Int sample_num_threshold  #minimum number of samples to get reliable joint calling results
  }
  command {
    set -euo pipefail
    if [[ $(cut -f 1 ~{sample_name_map} | wc -l) -ne $(cut -f 1 ~{sample_name_map} | sort | uniq | wc -l) ]]
    then
      echo "Samples in the sample_name_map are not unique" 1>&2
      exit 1
    elif [[ $(cut -f 1 ~{sample_name_map} | wc -l) -lt ~{sample_num_threshold} ]]
    then
      echo "There are fewer than ~{sample_num_threshold} samples in the sample_name_map" 1>&2
      echo "Having fewer than ~{sample_num_threshold} samples means there likely isn't enough data to complete joint calling" 1>&2
      exit 1
    else
      echo true
    fi
  }

  output {
    Boolean samples_unique = read_boolean(stdout())
  }

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "1000"
  }
}

############
### Split genome intervals for scattering
############
task SplitIntervalList {
  input {
    File GATK
    File interval_list
    Int scatter_count
    File ref_fasta
    File ref_index
    File ref_dict
    Boolean sample_names_unique_done
    String scatter_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
  }
  command {
    java -Xms3g \
      -jar ~{GATK} \
      SplitIntervals \
      -L ~{interval_list} -O scatterDir -scatter ~{scatter_count} -R ~{ref_fasta} \
      -imr OVERLAPPING_ONLY \
      -mode ~{scatter_mode}
   }

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "3500"
    runtime_minutes: "20"
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

############
### Import single-sample GVCFs into GenomicsDB before joint genotyping
############
task ImportGVCFs {
  input {
    File GATK
    File sample_name_map
    File interval
    File ref_fasta
    File ref_index
    File ref_dict
    String workspace_dir_name
    String TMP_DIR
    Int batch_size
  }
    # Using a nightly version of GATK containing fixes for GenomicsDB
    # https://github.com/broadinstitute/gatk/pull/5899

# /!\ if existing directory with same name, will be deleted

    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    # if export TILEDB_DISABLE_FILE_LOCKING=1 fails, need to upgrade to GATK4.1.8.0, which includes GenomicsDB (1.3.0)
    # 02.07 changes: deactivate consolidate option, make use of tmp_dir, increased memory
  command <<<
    set -euo pipefail

    rm -rf ~{workspace_dir_name} 
    TMP_DIR=`mktemp -d /tmp/tmp.XXXXXX`
    export TILEDB_DISABLE_FILE_LOCKING=1 #required when working on a POSIX filesystem (e.g. Lustre, NFS, xfs, ext4) before running any GenomicsDB tool
    cp ~{sample_name_map} tmp_map
    java -Xms16g -Djava.io.tmpdir=${TMP_DIR} \
      -jar ~{GATK} \
      GenomicsDBImport \
      --genomicsdb-workspace-path ~{workspace_dir_name} \
      --batch-size ~{batch_size} \
      -L ~{interval} \
      --sample-name-map tmp_map \
      --reader-threads 5 \
      --merge-input-intervals \
      --max-num-intervals-to-import-in-parallel 3 \
      --tmp-dir=${TMP_DIR}

    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
  >>>

  runtime {
    cpus: "4"
	  requested_memory_mb_per_core: "40000"
    runtime_minutes: "4500"
  }

  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}

############
### Perform joint genotyping on one or more samples
############
task GenotypeGVCFs {
  input {
    File GATK
    File workspace_tar
    File interval

    String output_vcf_filename

    File ref_fasta
    File ref_index
    File ref_dict
    File dbsnp_vcf # originally set as string???
    File dbsnp_vcf_index

    # This is needed for gVCFs generated with GATK3 HaplotypeCaller
    Boolean allow_old_rms_mapping_quality_annotation_data = false
  }
  command <<<

    set -euo pipefail

    tar -xf ~{workspace_tar}
    WORKSPACE=$(basename ~{workspace_tar} .tar)
    export TILEDB_DISABLE_FILE_LOCKING=1
    java -Xms16g \
      -jar ~{GATK} \
      GenotypeGVCFs \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      -D ~{dbsnp_vcf} \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --only-output-calls-starting-in-intervals \
      --use-new-qual-calculator \
      -V gendb://$WORKSPACE \
      -L ~{interval} \
      ~{true='--allow-old-rms-mapping-quality-annotation-data' false='' allow_old_rms_mapping_quality_annotation_data} \
      --merge-input-intervals
  >>>

  runtime {
    cpus: "2"
	  requested_memory_mb_per_core: "40000"  
    runtime_minutes: "5700"  
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

############
### Perform joint genotyping on one or more samples for large samplesets
############
task GnarlyGenotyper {
  input {
    File GATK
    File workspace_tar
    File interval
    String output_vcf_filename
    File ref_fasta
    File ref_index
    File ref_dict
    String dbsnp_vcf
  }
  command <<<
    set -e

    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)

    # use a query.json to set some params that aren't exposed -- ewwwww
    cat <<EOF > $WORKSPACE/query.json
      {
        "scan_full": true,
        "workspace": "genomicsdb",
        "array": "genomicsdb_array",
        "vid_mapping_file": "genomicsdb/vidmap.json",
        "callset_mapping_file": "genomicsdb/callset.json",
        "reference_genome": "/cromwell_root/broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
        "max_diploid_alt_alleles_that_can_be_genotyped": 6,
        "produce_GT_field": true
      }
    EOF

    java -Xms8g \
      -jar ~{GATK} \
      GnarlyGenotyper \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      --output-database-name annotationDB.vcf.gz \
      -D ~{dbsnp_vcf} \
      --only-output-calls-starting-in-intervals \
      --use-new-qual-calculator \
      -V gendb://$WORKSPACE \
      -L ~{interval} \
      -stand-call-conf 10 \
      --merge-input-intervals
  >>>

  runtime {
    cpus: "2"
	  requested_memory_mb_per_core: "26000"    
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
    File output_database = "annotationDB.vcf.gz"
    File output_database_index = "annotationDB.vcf.gz.tbi"
  }
}

############
### Filter variant calls using ExcessHet - EXCLUDE for samplesets <100s
############

## Hard-filter a large cohort callset on ExcessHet applies ONLY to callsets with a large number of samples, e.g. 100s. 
## Small cohorts should NOT trigger ExcessHet filtering as values should remain small. 
## This annotation estimates the probability of the called samples exhibiting excess heterozygosity with respect 
## to the null hypothesis that the samples are unrelated. The higher the score, the higher the chance that the variant 
## is a technical artifact or that there is consanguinuity among the samples.

task HardFilterAndMakeSitesOnlyVcf {
  input {
    File GATK
    File vcf
    File vcf_index
    Float excess_het_threshold

    String variant_filtered_vcf_filename
    String sites_only_vcf_filename
  }
#GATK? Annotated as picard tool in documentation
  command <<<
    set -euo pipefail

    java -Xms3g \
      -jar ~{GATK} \
      VariantFiltration \
      --filter-expression "ExcessHet > ~{excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ~{variant_filtered_vcf_filename} \
      -V ~{vcf}

    java -Xms3g \
      -jar ~{GATK} \
      MakeSitesOnlyVcf \
      -I ~{variant_filtered_vcf_filename} \
      -O ~{sites_only_vcf_filename}
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "4000"  
    runtime_minutes: "60" 
  }

  output {
    File variant_filtered_vcf = "~{variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "~{variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "~{sites_only_vcf_filename}"
    File sites_only_vcf_index = "~{sites_only_vcf_filename}.tbi"
  }
}

############
### Make sites Only - For samplesets <100s, when not hardfiltering
############

task MakeSitesOnlyVcf {
  input {
    File GATK
    File vcf
    File vcf_index
    String sites_only_vcf_filename
  }
#GATK? Annotated as picard tool in documentation
  command <<<
    set -euo pipefail
    java -Xms3g \
      -jar ~{GATK} \
      MakeSitesOnlyVcf \
      -I ~{vcf} \
      -O ~{sites_only_vcf_filename}
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "4000"   
    runtime_minutes: "60" 
  }

  output {
    File sites_only_vcf = "~{sites_only_vcf_filename}"
    File sites_only_vcf_index = "~{sites_only_vcf_filename}.tbi"
  }
}

############
### Gathers multiple VCF file
############
task GatherVcfs {
  input {
    File GATK
    File tabix
    Array[File] input_vcfs
    String output_vcf_name
  }
    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    # check that output combined vcf is sorted, else tabix will fail
  command <<<
    set -euo pipefail
    module add UHTS/Analysis/samtools/1.10
    java -Xms6g \
      -jar ~{GATK} \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ~{sep=" --input " input_vcfs} \
      --output tmp_vcf.gz
    
    bcftools sort tmp_vcf.gz  -O z -o ~{output_vcf_name}

    ~{tabix} ~{output_vcf_name}
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "9000" 
    runtime_minutes: "800"
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

############
### Indel recalibration model
############
task IndelsVariantRecalibrator {
  input {
    File GATK
    String recalibration_filename
    String tranches_filename

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File mills_resource_vcf
    File axiomPoly_resource_vcf
    File dbsnp_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 4
  }
## If a dataset gives <4 clusters, e.g. as can happen for smaller data, 
## then the tool will tell there is insufficient data with a No data found error message. 
## In this case, try decrementing the --max-gaussians value. 

  command <<<
    set -euo pipefail

    java -Xms24g \
      -jar ~{GATK} \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode INDEL \
      --max-gaussians ~{max_gaussians} \
      -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    cpus: "2"
	  requested_memory_mb_per_core: "26000"  
    runtime_minutes: "300"  
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}

############
### SNP recalibration model - parallel when >100k samples
############
task SNPsVariantRecalibratorCreateModel {
  input {
    File GATK
    String recalibration_filename
    String tranches_filename
    Int downsampleFactor
    String model_report_filename

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 6
  }

  command <<<
    set -euo pipefail

    java -Xms100g \
      -jar ~{GATK} \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode SNP \
      --sample-every-Nth-variant ~{downsampleFactor} \
      --output-model ~{model_report_filename} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    cpus: "2"
	  requested_memory_mb_per_core: "104000"   #104 Gib in original!
    runtime_minutes: "300"
  }

  output {
    File model_report = "~{model_report_filename}"
  }
}

############
### SNP recalibration model - standard
############
task SNPsVariantRecalibrator {
  input {
    File GATK
    String recalibration_filename
    String tranches_filename
    File? model_report
    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values
    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index
    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 6
    Int? machine_mem_gb
  }
  Int auto_mem = ceil(2 * size([sites_only_variant_filtered_vcf, hapmap_resource_vcf, omni_resource_vcf,one_thousand_genomes_resource_vcf, dbsnp_resource_vcf],"GiB"))
  Int machine_mem = select_first([machine_mem_gb, if auto_mem < 7 then 7 else auto_mem])
  Int java_mem = machine_mem - 1

  String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

  command <<<
    set -euo pipefail

    MODEL_REPORT=~{model_report}

    java -Xms~{java_mem}g \
      -jar ~{GATK} \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode SNP \
      ~{model_report_arg} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    cpus: "2"
	  requested_memory_mb_per_core: "~{machine_mem}000"  
    runtime_minutes: "300"
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}
## Each recal model step produces a .recal recalibration table and a .tranches tranches table. 
## In the filtering step, ApplyVQSR will use both types of data. 


############
### Filter variants using VQSR recal model
############
task ApplyRecalibration {
 input {
    File GATK
    String recalibrated_vcf_filename
    File input_vcf
    File input_vcf_index
    File indels_recalibration
    File indels_recalibration_index
    File indels_tranches
    File snps_recalibration
    File snps_recalibration_index
    File snps_tranches
    Float indel_filter_level
    Float snp_filter_level
    Boolean use_allele_specific_annotations
  }
  command <<<
    set -euo pipefail

    java -Xms5g \
      -jar ~{GATK} \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ~{input_vcf} \
      --recal-file ~{indels_recalibration} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      --tranches-file ~{indels_tranches} \
      --truth-sensitivity-filter-level ~{indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    java -Xms5g \
      -jar ~{GATK} \
      ApplyVQSR \
      -O ~{recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ~{snps_recalibration} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      --tranches-file ~{snps_tranches} \
      --truth-sensitivity-filter-level ~{snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "7000" 
    runtime_minutes: "200" 
  }

  output {
    File recalibrated_vcf = "~{recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "~{recalibrated_vcf_filename}.tbi"
  }
}

############
### Collect variant metrics
############
task CollectVariantCallingMetrics {
  input {
    File GATK
    File input_vcf
    File input_vcf_index
    String metrics_filename_prefix
    File dbsnp_vcf
    File dbsnp_vcf_index
    File interval_list
    File ref_dict
  }
#GATK? Picard tool in documentation
  command <<<
    set -euo pipefail

    java -Xms6g \
      -jar ~{GATK} \
      CollectVariantCallingMetrics \
      --INPUT ~{input_vcf} \
      --DBSNP ~{dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ~{ref_dict} \
      --OUTPUT ~{metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ~{interval_list}
  >>>

  output {
    File detail_metrics_file = "~{metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "~{metrics_filename_prefix}.variant_calling_summary_metrics"
  }

  runtime {
    cpus: "2"
	  requested_memory_mb_per_core: "8000" 
    runtime_minutes: "500" 
  }
}

############
### Select variants sites from Fingerprint
############
task SelectFingerprintSiteVariants {
  input {
    File GATK
    File input_vcf
    File haplotype_database
    String base_output_name
  }
  command <<<
    set -euo pipefail

    function hdb_to_interval_list() {
        input=$1
        awk 'BEGIN{IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {print $1,$2,$2,"+","interval-"NR}' $1
    }

    hdb_to_interval_list ~{haplotype_database} > hdb.interval_list

    java -Xms6g \
      -jar ~{GATK} \
      SelectVariants \
      --variant ~{input_vcf} \
      --intervals hdb.interval_list \
      --output ~{base_output_name}.vcf.gz
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "8000" 
  }

  output {
    File output_vcf = "~{base_output_name}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.vcf.gz.tbi"
  }
}

############
### Cross Check fingerprint
############
task CrossCheckFingerprint {
  input {
    File PICARD
    Array[File] gvcf_paths
    Array[File] gvcf_index_paths
    Array[File] vcf_paths
    Array[File] vcf_index_paths
    File sample_name_map
    File haplotype_database
    String output_base_name
    Boolean scattered = false
    Array[String] expected_inconclusive_samples = []
  }
  Int num_gvcfs = length(gvcf_paths)
  Int cpu = if num_gvcfs < 32 then num_gvcfs else 32
  # Compute memory to use based on the CPU count, following the pattern of
  # 3.75GiB / cpu used by GCP's pricing: https://cloud.google.com/compute/pricing
  Int memMb = round(cpu * 3.75 * 1024)
  String output_name = output_base_name + ".fingerprintcheck"

#      --SKIP_INPUT_READABLITY_TEST \  unrecognized option?
  command <<<
    set -eu

    gvcfInputsList=~{write_lines(gvcf_paths)}
    vcfInputsList=~{write_lines(vcf_paths)}

    cp $gvcfInputsList gvcf_inputs.list
    cp $vcfInputsList vcf_inputs.list

    java -Dpicard.useLegacyParser=false -Xms~{memMb - 512}m \
      -jar ~{PICARD} \
      CrosscheckFingerprints \
      --INPUT gvcf_inputs.list \
      --SECOND_INPUT vcf_inputs.list \
      --HAPLOTYPE_MAP ~{haplotype_database} \
      --INPUT_SAMPLE_FILE_MAP ~{sample_name_map} \
      --CROSSCHECK_BY SAMPLE \
      --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
      --NUM_THREADS ~{cpu} \
      ~{true='--EXIT_CODE_WHEN_MISMATCH 0' false='' scattered} \
      --OUTPUT ~{output_name}
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "~{memMb}" 
  }

  output {
    File crosscheck_metrics = output_name
  }
}

############
### Gather metrics
############
task GatherPicardMetrics {
  input {
    Array[File] metrics_files
    String output_file_name
  }
  # Don't use this task to gather tens of thousands of files.
  # Cromwell can't handle it.
  # This cannot gather metrics with histograms
  command {
    head -n 7 ~{metrics_files[0]} > ~{output_file_name}
    for metrics_file in ~{sep=' ' metrics_files}; do
      sed -n '1,7d;p' $metrics_file | grep -v '^$' >> ~{output_file_name}
    done
  }

  output {
    File gathered_metrics = "~{output_file_name}"
  }

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "4000" 
  }
}

############
### Refines genotype quality using all cohort samples (if >10) and the 1000G dataset
############# 
task CalculatGenotypePosteriors {
  input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File vcf
    File vcf_index
    File one_thousand_genomes_vcf
    File one_thousand_genomes_vcf_index
    String output_vcf_filename
  }
  command {
    java -Xmx4g \
    -jar ~{GATK} \
    CalculateGenotypePosteriors \
    -R ~{ref_fasta} \
    --supporting ~{one_thousand_genomes_vcf} \
    -V ~{vcf} \
    -O ~{output_vcf_filename}
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "6000"
  runtime_minutes: "1000"
	queue: "normal"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

############
### Annotates variants with refined GQ <20 and DP<10 as lowQUAL
#############  

task VariantFilterLowQ {
 input {
    File GATK
    File ref_fasta
    File ref_index
    File ref_dict
    File vcf
    File vcf_index
    String output_vcf_filename
  }
  command {
    java -Xms3g \
    -jar ~{GATK} \
    VariantFiltration \
    -R ~{ref_fasta} \
    -V ~{vcf} \
    -G-filter "GQ < 20.0" \
    -G-filter-name lowGQ \
    -G-filter "DP < 10.0" \
    -G-filter-name lowDP \
    -O ~{output_vcf_filename}
  }
  runtime {
	cpus: "1"
	requested_memory_mb_per_core: "4000"
  runtime_minutes: "1000"
	queue: "normal"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

############
### Create sampleset partitions
############
task PartitionSampleNameMap {
  input {
    File sample_name_map
    Int line_limit
  }
# Let the OS catch up with creation of files for glob command
  command {
    cut -f 2 ~{sample_name_map} > sample_paths
    split -l ~{line_limit} -d sample_paths partition_
    sleep 1
  }

  output {
    Array[File] partitions = glob("partition_*")
  }

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "1000" 
  }
}

############
### Annotate vcf with annovar (small contigs)
############
task AnnovarScatteredVCF {
  input {
    File input_vcf
    String AnnovarDB   ##path to AnnovarDB (hg19/hg38)
    String genome_build ## hg38 / hg18 / hg19. Needs corresponding dbs to be downloaded!
    String base_vcf_name  
  }

  command <<<
  module add UHTS/Analysis/EPACTS/3.2.6
  export PATH="$PATH:/data/PRTNR/CHUV/MED/jfellay/default_sensitive/redin/tools/annovar"
  table_annovar.pl ~{input_vcf} \
  "~{AnnovarDB}" -buildver ~{genome_build} \
  -out ~{base_vcf_name} -remove \
  -protocol refGene,ensGene,gnomad30_genome,dbnsfp41a,1000g2015aug_all,kaviar_20150923,clinvar,cytoBand,dbscsnv11,spidex_lifted \
  -operation g,g,f,f,f,f,f,r,f,f -nastring . -vcfinput --thread 2 --polish

  bgzip "~{base_vcf_name}.hg38_multianno.vcf"
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "14000"
    runtime_minutes: "1500"
  }

  output {
    File output_vcf = "~{base_vcf_name}.hg38_multianno.vcf.gz"
  }
}

############
### Annotate vcf with vcfanno (small contigs)
############
task vcfannoScatteredVCF {
  input {
    File input_vcf
    File conf_file
    String base_vcf_name  
  }

  command <<<
  module add UHTS/Analysis/vcfanno/0.3.2
  module add UHTS/Analysis/EPACTS/3.2.6
  vcfanno ~{conf_file} ~{input_vcf} > "~{base_vcf_name}.hg38_vcfanno.vcf"

  bgzip "~{base_vcf_name}.hg38_vcfanno.vcf"
  tabix "~{base_vcf_name}.hg38_vcfanno.vcf.gz"
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "6000"
    runtime_minutes: "60"
  }

  output {
    File output_vcf = "~{base_vcf_name}.hg38_vcfanno.vcf.gz"
    File output_vcf_index = "~{base_vcf_name}.hg38_vcfanno.vcf.gz.tbi"
  }
}

############
### Annotate vcf with VEP
############
task VEPannoScatteredVCF {
  input {
    File input_vcf
    String base_vcf_name  
  }

 #pick option: to select annotations from canonical transcripts
  command <<<
  source /dcsrsoft/spack/bin/setup_dcsrsoft
  module load gcc
  module load singularity
  module load htslib/1.12

  export PERL5LIB=$PERL5LIB:/db/local/vep/Plugins.  # pour spécifier le répertoire des plugins
  export SINGULARITY_BINDPATH="/scratch,/users,/dcsrsoft,/db"

  singularity run /dcsrsoft/singularity/containers/ensembl-vep_104.sif vep \
  -i ~{input_vcf} \
  --plugin dbNSFP,/db/local/vep/plugins_data/dbNSFP4.1a_grch38.gz,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
  --plugin SpliceAI,snv=/db/local/vep/plugins_data/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/db/local/vep/plugins_data/spliceai_scores.raw.indel.hg38.vcf.gz \
  --buffer_size 100000 \
  --offline --fork 10 \
  --dir_cache=/db/local \
  --vcf --force_overwrite \
  --pick --cache -o "~{base_vcf_name}.hg38.VEP.vcf"

  bgzip "~{base_vcf_name}.hg38.VEP.vcf"
  tabix "~{base_vcf_name}.hg38.VEP.vcf.gz"
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "20000"
    runtime_minutes: "2000"
  }

  output {
    File output_vcf = "~{base_vcf_name}.hg38.VEP.vcf.gz"
    File output_vcf_index = "~{base_vcf_name}.hg38.VEP.vcf.gz.tbi"
  }
}

############
### Split vcf
############
task SplitVCF {
  input {
    File input_vcf
    File input_vcf_index
    File interval
    String base_vcf_name  
  }

  command <<<
  module add UHTS/Analysis/samtools/1.10
  module add UHTS/Analysis/EPACTS/3.2.6
  grep '^chr' ~{interval} > "~{interval}_clean"
  bcftools view -R "~{interval}_clean" ~{input_vcf} -O z -o "~{base_vcf_name}.hg38.vcf.gz"
  tabix "~{base_vcf_name}.hg38.vcf.gz"
  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "6000"
    runtime_minutes: "200"
  }

  output {
    File output_vcf = "~{base_vcf_name}.hg38.vcf.gz"
  }
}

############
### Annotate vcf with spliceAI
############
task spliceAIScatteredVCF {
  input {
    File input_vcf
    File conf_file
    String base_vcf_name  
  }

  command <<<
  spliceai -I input.vcf -O output.vcf -R genome.fa -A grch37
  cat input.vcf | spliceai -R genome.fa -A grch37 > output.vcf

  >>>

  runtime {
    cpus: "1"
	  requested_memory_mb_per_core: "6000"
    runtime_minutes: "20"
  }

  output {
    File output_vcf = "~{base_vcf_name}.hg38_vcfanno.vcf.gz"
  }
}