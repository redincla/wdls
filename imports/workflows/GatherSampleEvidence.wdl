##########################################################################################

## Base script:   https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherSampleEvidence.wdl

##########################################################################################

version 1.0

import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/tasks/SVtasks.wdl" as SVtasks
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/workflows/runDelly-GATKSV.wdl" as delly
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/workflows/runManta-GATKSV.wdl" as manta
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/workflows/runSmoove-GATKSV.wdl" as smoove
import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/workflows/CramToBam.ReviseBase.wdl" as ctb_revise
#import "/scratch/beegfs/PRTNR/CHUV/MED/jfellay/default_sensitive/WGS/wdls/imports/workflows/rubWhamg.wdl" as wham

workflow GatherSampleEvidence {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id

    # Use to revise Y, R, W, S, K, M, D, H, V, B, X bases in BAM to N. Use only if providing a CRAM file as input 
    Boolean? revise_base_cram_to_bam
    File? primary_contigs_fai # required if using revise_base_cram_to_bam (or if run_module_metrics = true)

    # Evidence collection flags
    Boolean collect_coverage = true
    Boolean collect_pesr = true

    # Common parameters
    File primary_contigs_list
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta
    String? reference_version   # Either "38" or "19"

    # Coverage collection inputs
    File preprocessed_intervals

    # Delly specific inputs
    Boolean? run_delly 
    File Delly_exclude_regions_bed

    # Manta specific inputs
    Boolean? run_manta
    File Manta_region_bed
    File Manta_region_bed_index

    # smoove inputs
    Boolean? run_smoove
    File Smoove_exclude_regions_bed
  }
    
    Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)

  # Convert to BAM if we have a CRAM
  if (!is_bam_) {
        call SVtasks.FixHeaderCRAM {  #fix time field in cram header for later smoove run
            input:
                input_cram = bam_or_cram_file,
                input_cram_index = bam_or_cram_index
        }

        call ctb_revise.CramToBamReviseBase {
            input:
                input_cram = FixHeaderCRAM.output_cram,
                input_cram_index = FixHeaderCRAM.output_cram_index,
                ref_fasta = reference_fasta,
                ref_fasta_fai = reference_index,
                contiglist = select_first([primary_contigs_fai]),
                base_name = sample_id
        }
    }

    if (is_bam_) {
        call SVtasks.FixHeaderBAM {  #fix time field in cram header for later smoove run
            input:
                input_bam = bam_or_cram_file,
                input_bam_index = bam_or_cram_index
        }
    }

  File bam_file = select_first([CramToBamReviseBase.bam_file, FixHeaderBAM.output_bam])
  File bam_index = select_first([CramToBamReviseBase.bam_index, FixHeaderBAM.output_bam_index])

   call SVtasks.CollectCounts {
      input:
        intervals = preprocessed_intervals,
        bam = bam_file,
        bam_idx = bam_index,
        sample_id = sample_id,
        ref_fasta = reference_fasta,
        ref_fasta_fai = reference_index,
        ref_fasta_dict = reference_dict
    }

    call delly.Delly {
      input:
        bam_or_cram_file = bam_file,
        bam_or_cram_index = bam_index,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        exclude_regions_bed = Delly_exclude_regions_bed
    }

    call manta.Manta {
      input:
        bam_or_cram_file = bam_file,
        bam_or_cram_index = bam_index,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        region_bed = Manta_region_bed,
        region_bed_index = Manta_region_bed_index
    }

    call smoove.Smoove {
      input:
        bam_or_cram_file = bam_file,
        bam_or_cram_index = bam_index,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        region_bed = Smoove_exclude_regions_bed
    }

    call SVtasks.PESRCollection {
      input:
        cram = bam_file,
        cram_index = bam_index,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        reference_dict = reference_dict
    }

  output {
    File coverage_counts = CollectCounts.counts

    File delly_vcf = Delly.vcf
    File delly_index = Delly.index

    File manta_vcf = Manta.vcf
    File manta_index = Manta.index

    File smoove_vcf = Smoove.vcf
    File smoove_index = Smoove.index

    File pesr_disc = PESRCollection.disc_out
    File pesr_disc_index = PESRCollection.disc_out_index
    File pesr_split = PESRCollection.split_out
    File pesr_split_index = PESRCollection.split_out_index
  }
}