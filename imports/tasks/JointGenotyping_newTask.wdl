############
### Add new samples to existing GenomicsDB before joint genotyping
############
task updateGenomicsDB {
  input {
    File GATK
    File sample_name_map
    File workspace_tar
    String TMP_DIR
    String output_interval_file_name
    Int batch_size
  }
    # The --genomicsdb-update-workspace-path must point to a existing genomicsdb workspace
    # /!\ Always backup existing genomicsdb workspaces before adding new samples. 
    # If the tool fails during incremental import for any reason, the workspace may be corrupted

  command <<<
    set -euo pipefail

    TMP_DIR=`mktemp -d /tmp/tmp.XXXXXX`
    export TILEDB_DISABLE_FILE_LOCKING=1 #required when working on a POSIX filesystem (e.g. Lustre, NFS, xfs, ext4) before running any GenomicsDB tool
    cp ~{sample_name_map} tmp_map
    #need untar'd db
    tar -xf ~{workspace_tar}
    WORKSPACE=$(basename ~{workspace_tar} .tar)
    java -Xms16g -Djava.io.tmpdir=${TMP_DIR} \
      -jar ~{GATK} \
      GenomicsDBImport \
      --genomicsdb-update-workspace-path ${WORKSPACE} \
      --batch-size ~{batch_size} \
      --sample-name-map tmp_map \
      --reader-threads 5 \
      --tmp-dir=${TMP_DIR}

    #cannnot export interval list and update db at the same time
    java -Xms16g -Djava.io.tmpdir=${TMP_DIR} \
      -jar ~{GATK} \
      GenomicsDBImport \
      --genomicsdb-update-workspace-path ${WORKSPACE} \
      --output-interval-list-to-file ~{output_interval_file_name} \

    tar -cf ${WORKSPACE}.tar ${WORKSPACE}
  >>>

  runtime {
    cpus: "4"
	  requested_memory_mb_per_core: "40000"
    runtime_minutes: "2500"
  }

  output {
    File output_genomicsdb = "${WORKSPACE}.tar"
    File interval = "~{output_interval_file_name}"
  }
}

