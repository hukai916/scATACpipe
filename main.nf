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
def checkPathParamList = [ params.input_fragment, params.input_fastq, params.ref_bwa_index, params.ref_minimap2_index, params.ref_cellranger_index, params.ref_gtf, params.ref_fasta, params.whitelist_barcode, params.archr_genome_fasta, params.archr_blacklist, params.archr_scrnaseq ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters: already in initialise
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////

// Parse samplesheet:
if (params.input_fastq) {
  Channel
  .from(file(params.input_fastq, checkIfExists: true))
  .splitCsv(header: true, sep: ",", strip: true)
  .map {
    row ->
      [ row.sample_name, row.path_fastq_1, row.path_fastq_2, row.path_barcode ]
  }
  .unique()
  .set { ch_samplesheet_preprocess }
} else if (params.input_fragment) { // Parse ArchR samplesheet:
  Channel
  .from(file(params.input_fragment, checkIfExists: true))
  .splitCsv(header: true, sep: ",", strip: true)
  .map {
    row ->
      [ row.sample_name, row.file_path ]
  }
  .unique()
  .set { ch_samplesheet_archr }
} else {
  exit 1, "Must specify eitehr --input_fragment or --input_fastq!"
}

// Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////
include { PREPROCESS_DEFAULT } from './workflows/preprocess_default'
include { PREPROCESS_10XGENOMICS } from './workflows/preprocess_10xgenomics'
include { CHROMAP_INDEX } from './modules/local/chromap_index' addParams( options: modules['chromap_index'] )
include { CHROMAP_ATAC } from './modules/local/chromap_atac' addParams( options: modules['chromap_atac'] )
include { GET_WHITELIST_CHROMAP } from './modules/local/get_whitelist_chromap' addParams( options: modules['get_whitelist_chromap'] )

// include { PREPROCESS_CHROMAP } from './workflows/preprocess_chromap'
include { DOWNSTREAM_ARCHR } from './workflows/downstream_archr'
include { SPLIT_BED  } from './modules/local/split_bed' addParams( options: modules['split_bed'] )
include { SPLIT_BAM  } from './modules/local/split_bam' addParams( options: modules['split_bam'] )
include { MULTIQC    } from './modules/local/multiqc' addParams( options: modules['multiqc'] )

include { INPUT_CHECK_FRAGMENT } from './subworkflows/local/input_check_fragment'
include { INPUT_CHECK_FASTQ } from './subworkflows/local/input_check_fastq'

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)

workflow SCATACPIPE {
  take:
    input_fragment
    input_fastq

  main:
    if (input_fragment) {
      log.info "Running DOWNSTREAM ..."
      log.info "Validating sample sheet ... If pipeline exits, check .nextflow.log file."
      INPUT_CHECK_FRAGMENT (Channel.fromPath(input_fragment))

      DOWNSTREAM_ARCHR (INPUT_CHECK_FRAGMENT.out.fragment, "preprocess_null", "token1", "token2", "token3", "token4", "token5", "token6")
      SPLIT_BED (DOWNSTREAM_ARCHR.out[1])
      MULTIQC (DOWNSTREAM_ARCHR.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
    } else if (input_fastq) {
      log.info "Running PREPROCESS + DOWNSTREAM ..."
      log.info "Validating sample sheet ... If pipeline exits, check .nextflow.log file."
      INPUT_CHECK_FASTQ (Channel.fromPath(input_fastq))

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
        PREPROCESS_DEFAULT (INPUT_CHECK_FASTQ.out.reads, INPUT_CHECK_FASTQ.out.sample_count)
        DOWNSTREAM_ARCHR (PREPROCESS_DEFAULT.out[2], "preprocess_default", prep_genome_run, PREPROCESS_DEFAULT.out[5], PREPROCESS_DEFAULT.out[6], prep_gtf_run, PREPROCESS_DEFAULT.out[7], PREPROCESS_DEFAULT.out[8])
        SPLIT_BED (DOWNSTREAM_ARCHR.out[1]) // take a tuple (sample_name, fragment_path, tsv_path) as input
        SPLIT_BAM (PREPROCESS_DEFAULT.out[3], DOWNSTREAM_ARCHR.out[2].collect(), PREPROCESS_DEFAULT.out[4].collect(), params.barcode_regex) // input: sample_name, all_bams, all_fragments, barcode_regex
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
        PREPROCESS_10XGENOMICS (INPUT_CHECK_FASTQ.out.reads, INPUT_CHECK_FASTQ.out.sample_count)
        // DOWNSTREAM_ARCHR (PREPROCESS_10XGENOMICS.out[2], "preprocess_10xgenomics")
        DOWNSTREAM_ARCHR (PREPROCESS_10XGENOMICS.out[2], "preprocess_10xgenomics", prep_genome_run, PREPROCESS_10XGENOMICS.out[5], PREPROCESS_10XGENOMICS.out[6], prep_gtf_run, PREPROCESS_10XGENOMICS.out[7], PREPROCESS_10XGENOMICS.out[8])
        SPLIT_BED (DOWNSTREAM_ARCHR.out[1])
        SPLIT_BAM (PREPROCESS_10XGENOMICS.out[3], DOWNSTREAM_ARCHR.out[2].collect(), PREPROCESS_10XGENOMICS.out[4].collect(), "NA")
        MULTIQC (DOWNSTREAM_ARCHR.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
      } else if (params.preprocess == "chromap") {
        // Determine if PREP_GENOME and PREP_GTF run or not_run:
        // // If index folder supplied: both PREP_GENOME and PREP_GTF must not_run
        prep_genome_run = "run"
        prep_gtf_run    = "run"
        if (params.ref_chromap_index) {
          prep_genome_run = "not_run"
          prep_gtf_run    = "not_run"
        }

        // PREPROCESS_10XGENOMICS (ch_samplesheet)
        PREPROCESS_CHROMAP (INPUT_CHECK_FASTQ.out.reads, INPUT_CHECK_FASTQ.out.sample_count)
        // DOWNSTREAM_ARCHR (PREPROCESS_10XGENOMICS.out[2], "preprocess_10xgenomics")
        DOWNSTREAM_ARCHR (PREPROCESS_CHROMAP.out[2], "preprocess_10xgenomics", prep_genome_run, PREPROCESS_10XGENOMICS.out[5], PREPROCESS_10XGENOMICS.out[6], prep_gtf_run, PREPROCESS_10XGENOMICS.out[7], PREPROCESS_10XGENOMICS.out[8])
        SPLIT_BED (DOWNSTREAM_ARCHR.out[1])
        // SPLIT_BAM (PREPROCESS_10XGENOMICS.out[3], DOWNSTREAM_ARCHR.out[2].collect(), PREPROCESS_10XGENOMICS.out[4].collect(), "NA")
        MULTIQC (DOWNSTREAM_ARCHR.out[0].ifEmpty([]).mix(Channel.from(ch_multiqc_config)).collect())
      } else {
        exit 1, "must supply valid --preprocess option"
      }
    } else {
      exit 1, "Pls supply either --input_fragment or --input_fastq"
    }
}

workflow {
  // SCATACPIPE (params.input_fragment, params.input_fastq, ch_samplesheet_preprocess)
  SCATACPIPE (params.input_fragment, params.input_fastq)


  // if (params.input_fragment) {
  //   // SCATACPIPE (params.input_fragment, params.input_fastq, ch_samplesheet_archr)
  //   SCATACPIPE (params.input_fragment, params.input_fastq, params.input_fragment)
  // } else if (params.input_fastq) {
  //   SCATACPIPE (params.input_fragment, params.input_fastq, ch_samplesheet_preprocess)
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
