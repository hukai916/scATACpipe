#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/scatacpipe
========================================================================================
    Github : https://github.com/nf-core/scatacpipe
    Website: https://nf-co.re/scatacpipe
    Slack  : https://nfcore.slack.com/channels/scatacpipe
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { SCATACPIPE } from './workflows/scatacpipe'

//
// WORKFLOW: Run main nf-core/scatacpipe analysis pipeline
//
workflow NFCORE_SCATACPIPE {
    SCATACPIPE ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_SCATACPIPE ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
