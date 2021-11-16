#!/usr/bin/env nextflow
/*
========================================================================================
    scatacpipe
========================================================================================
    Github : https://github.com/hukai916/scatacpipe
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def modules = params.modules.clone()


/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

// Check input path parameters to see if they exist:
def checkPathParamList = [ params.input_archr, params.input_preprocess, params.ref_bwa_index, params.ref_minimap2_index, params.ref_cellranger_index, params.ref_gtf, params.ref_fasta, params.barcode_whitelist, params.archr_genome_fasta, params.archr_blacklist_bed, params.archr_scrnaseq ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters: already in initialise
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////

// Parse samplesheet:
if (params.input_preprocess) {
  Channel
  .from(file(params.input_preprocess, checkIfExists: true))
  .splitCsv(header: true, sep: ",", strip: true)
  .map {
    row ->
      [ row.sample_name, row.path_fastq_1, row.path_fastq_2, row.path_barcode ]
  }
  .unique()
  .set { ch_samplesheet_preprocess }
} else if (params.input_archr) { // Parse ArchR samplesheet:
  Channel
  .from(file(params.input_archr, checkIfExists: true))
  .splitCsv(header: true, sep: ",", strip: true)
  .map {
    row ->
      [ row.sample_name, row.file_path ]
  }
  .unique()
  .set { ch_samplesheet_archr }
} else {
  exit 1, "Must specify eitehr --input_archr or --input_preprocess!"
}

// Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////
include { PREPROCESS_DEFAULT } from './workflows/preprocess_default'
include { PREPROCESS_10XGENOMICS } from './workflows/preprocess_10xgenomics'
include { DOWNSTREAM_ARCHR } from './workflows/downstream_archr'
include { SPLIT_BED  } from './modules/local/split_bed' addParams( options: modules['split_bed'] )
include { SPLIT_BAM  } from './modules/local/split_bam' addParams( options: modules['split_bam'] )
include { MULTIQC    } from './modules/local/multiqc' addParams( options: modules['multiqc'] )

include { INPUT_CHECK_ARCHR } from './subworkflows/local/input_check_archr'
include { INPUT_CHECK_PREPROCESS } from './subworkflows/local/input_check_preprocess'

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)

workflow SCATACPIPE {
  take:
    input_archr
    input_preprocess
    // ch_samplesheet

  main:
    if (input_archr) {
      log.info "Running DOWNSTREAM ..."
      log.info "Validating sample sheet ... If pipeline exits, check .nextflow.log file."
      INPUT_CHECK_ARCHR (Channel.fromPath(input_archr))

      DOWNSTREAM_ARCHR (INPUT_CHECK_ARCHR.out.fragment, "preprocess_null", "token1", "token2", "token3", "token4", "token5", "token6")
      SPLIT_BED (DOWNSTREAM_ARCHR.out[1])
      MULTIQC (DOWNSTREAM_ARCHR.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
    } else if (input_preprocess) {
      log.info "Running PREPROCESS + DOWNSTREAM ..."
      log.info "Validating sample sheet ... If pipeline exits, check .nextflow.log file."
      INPUT_CHECK_PREPROCESS (Channel.fromPath(input_preprocess))

      if (params.preprocess == "default") {
        // Determine if PREP_GENOME and PREP_GTF run or not_run:
        //// PREP_GTF is always not_run
        //// If index folder supplied: PREP_GENOME must not run
        prep_genome_run = "run"
        prep_gtf_run    = "not_run"
        if ((params.mapper == "bwa") && params.ref_bwa_index) {
          prep_genome_run = "not_run"
        } else if ((params.mapper == "minimap2") && params.ref_minimap2_index) {
          prep_genome_run = "not_run"
        }

        // PREPROCESS_DEFAULT (ch_samplesheet)
        PREPROCESS_DEFAULT (INPUT_CHECK_PREPROCESS.out.reads)
        DOWNSTREAM_ARCHR (PREPROCESS_DEFAULT.out[2], "preprocess_default", prep_genome_run, PREPROCESS_DEFAULT.out[5], PREPROCESS_DEFAULT.out[6], prep_gtf_run, PREPROCESS_DEFAULT.out[7], PREPROCESS_DEFAULT.out[8])
        SPLIT_BED (DOWNSTREAM_ARCHR.out[1]) // take a tuple (sample_name, fragment_path, tsv_path) as input
        SPLIT_BAM (PREPROCESS_DEFAULT.out[3], DOWNSTREAM_ARCHR.out[2].collect(), PREPROCESS_DEFAULT.out[4].collect(), "[^:]*") // input: sample_name, all_bams, all_fragments, barcode_regex
        MULTIQC(PREPROCESS_DEFAULT.out[0].mix(DOWNSTREAM_ARCHR.out[0].ifEmpty([])).mix(Channel.from(ch_multiqc_config)).collect())
      } else if (params.preprocess == "10xgenomics") {
        // Determine if PREP_GENOME and PREP_GTF run or not_run:
        // // If index folder supplied: both PREP_GENOME and PREP_GTF must not_run
        prep_genome_run = "run"
        prep_gtf_run    = "run"
        if (params.ref_cellranger_index) {
          prep_genome_run = "not_run"
          prep_gtf_run    = "not_run"
        }

        // PREPROCESS_10XGENOMICS (ch_samplesheet)
        PREPROCESS_10XGENOMICS (INPUT_CHECK_PREPROCESS.out.reads)
        // DOWNSTREAM_ARCHR (PREPROCESS_10XGENOMICS.out[2], "preprocess_10xgenomics")
        DOWNSTREAM_ARCHR (PREPROCESS_10XGENOMICS.out[2], "preprocess_10xgenomics", prep_genome_run, PREPROCESS_10XGENOMICS.out[5], PREPROCESS_10XGENOMICS.out[6], prep_gtf_run, PREPROCESS_10XGENOMICS.out[7], PREPROCESS_10XGENOMICS.out[8])
        SPLIT_BED (DOWNSTREAM_ARCHR.out[1])
        SPLIT_BAM (PREPROCESS_10XGENOMICS.out[3], DOWNSTREAM_ARCHR.out[2].collect(), PREPROCESS_10XGENOMICS.out[4].collect(), "NA")
        MULTIQC (DOWNSTREAM_ARCHR.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
      } else if (params.preprocess == "biorad") {
        exit 1, "biorad to be added"
      } else {
        exit 1, "must supply valid --preprocess option"
      }
    } else {
      exit 1, "Pls supply either --input_archr or --input_preprocess"
    }
}

workflow {
  // SCATACPIPE (params.input_archr, params.input_preprocess, ch_samplesheet_preprocess)
  SCATACPIPE (params.input_archr, params.input_preprocess)


  // if (params.input_archr) {
  //   // SCATACPIPE (params.input_archr, params.input_preprocess, ch_samplesheet_archr)
  //   SCATACPIPE (params.input_archr, params.input_preprocess, params.input_archr)
  // } else if (params.input_preprocess) {
  //   SCATACPIPE (params.input_archr, params.input_preprocess, ch_samplesheet_preprocess)
  // }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
