version 1.0

## Copyright CHUV, 2021
## Script to launch gene-burden tests
## Following workflow from https://github.com/DrGBL/WES.WGS

#################################################################
# WORKFLOW DEFINITION - Gene Burden test of small case cohorts vs gnomAD
#################################################################

## Local import
import "./imports/workflows/BurdenTests.wdl" as CasesPrep
import "./imports/workflows/burdenTest-gnomAD.wdl" as CtrlsPrep
import "/home/credin/scratch/WGS/wdls/imports/tasks/JointCalling-tasks-test.wdl" as Tasks


#################################################################
# WORKFLOW DEFINITION
#################################################################

workflow BurdenTest {

String pipeline_version = "1.0"

input {
File cases_input_vcf
File cases_input_vcf_index
File ctrls_input_vcf
File ctrls_input_vcf_index
File ref_fasta
File ref_index
File ref_dict
File ped_file
File ancestry_IDs
String genome_version
File GATK
File tabix
File unpadded_intervals_file
Int scatter_count
Float AF_threshold
Int AC_threshold
File regions_list_cases
File regions_list_controls
}

String vcf_prefix=basename(cases_input_vcf, ".vcf.gz")

call CasesPrep.BurdenTest {
    input:
    input_vcf = cases_input_vcf,
    input_vcf_index = cases_input_vcf_index,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    ref_dict = ref_dict,
    ped_file = ped_file,
    ancestry_IDs = ancestry_IDs,
    genome_version = genome_version,
    GATK = GATK,
    tabix = tabix,
    unpadded_intervals_file = unpadded_intervals_file,
    scatter_count = scatter_count,
    regions_list = regions_list_cases,
    AF_threshold = AF_threshold,
    AC_threshold = AC_threshold
}

call CtrlsPrep.BurdenTestgnomAD {
    input:
    input_vcf = ctrls_input_vcf,
    input_vcf_index = ctrls_input_vcf_index,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    ref_dict = ref_dict,
    genome_version = genome_version,
    GATK = GATK,
    tabix = tabix,
    unpadded_intervals_file = unpadded_intervals_file,
    regions_list = regions_list_controls,
    scatter_count = scatter_count,
    AF_threshold = AF_threshold,
}

call CasesPrep.RunBurdenTest as RunBurdenTestHighImpact {
    input:
    cases_count_file = CasesPrep.HighImpactPASS_count,
    ctrls_count_file = CtrlsPrep.HighImpactPASS_count,
    cases_size = cases_size,
    ctrls_size = ctrls_size,
    base_output_name = vcf_prefix
    }

call CasesPrep.RunBurdenTest as RunBurdenTestModImpactNonMis {
    input:
    cases_count_file = CasesPrep.ModImpactNonMisPASS_count,
    ctrls_count_file = CtrlsPrep.ModImpactNonMisPASS_count,
    cases_size = cases_size,
    ctrls_size = ctrls_size,
    base_output_name = vcf_prefix
    }

call CasesPrep.RunBurdenTest as RunBurdenTestModImpactMissense {
    input:
    cases_count_file = CasesPrep.ModImpactMissensePASS_count,
    ctrls_count_file = CtrlsPrep.ModImpactMissensePASS_count,
    cases_size = cases_size,
    ctrls_size = ctrls_size,
    base_output_name = vcf_prefix
    }

call CasesPrep.RunBurdenTest as RunBurdenTestLowImpactSyn {
    input:
    cases_count_file = CasesPrep.LowImpactSynPASS_count,
    ctrls_count_file = CtrlsPrep.LowImpactSynPASS_count,
    cases_size = cases_size,
    ctrls_size = ctrls_size,
    base_output_name = vcf_prefix
    }

output {
    File output_burden_HighImpact = RunBurdenTestHighImpact.burden_file
    File output_burden_ModImpactNonMis = RunBurdenTestModImpactNonMis.burden_file
    File output_burden_ModImpactMissense = RunBurdenTestModImpactMissense.burden_file
    File output_burden_LowImpactSyn = RunBurdenTestLowImpactSyn.burden_file
}