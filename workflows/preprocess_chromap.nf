/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
// NfcoreSchema pops compilation error.

// Validate input parameters
// WorkflowScatacpipe.initialise(params, log)

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Modules: local
include { CHROMAP_INDEX } from '../modules/local/chromap_index' addParams( options: modules['chromap_index'] )
include { CHROMAP_ATAC } from '../modules/local/chromap_atac' addParams( options: modules['chromap_atac'] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow PREPROCESS_CHROMAP {
  take:
    reads
    sample_count

  main:
    // Examine if all required parameters supplied:
    if (!params.ref_chromap_index && !(params.ref_fasta && params.ref_gtf) && !params.ref_fasta_ensembl && !params.ref_fasta_ucsc) {
      msg = "Must supply one from below:\n" + "Option1:\n  --ref_chromap_index\n" + "Option2:\n  --ref_fasta\n  --ref_gtf\n" + "Option3:\n  --ref_fasta_ensembl\n" + "Option4:\n  --ref_fasta_ucsc\n"
      log.error msg
      exit 1, "EXIT!"
    }
    // Above is redundant to WorkflowMain::initialise()
    // module: staging sample_name in case of inconsistency in sample names
    MATCH_SAMPLE_NAME (reads, sample_count)

    if (params.ref_chromap_index) {
      // if cellranger index folder provided:
      log.info "Parameter --ref_chromap_index supplied, will use it as index file."
      CHROMAP_ATAC (MATCH_SAMPLE_NAME.out.sample_name.unique(), MATCH_SAMPLE_NAME.out.sample_files.collect(), params.ref_chromap_index)
    } else {
      if (params.ref_fasta) {
        if (params.ref_gtf) {
          log.info "Parameter --ref_fasta/ref_gtf supplied, will build index."
          // Module: prep_genome
          PREP_GENOME (params.ref_fasta, "custom_genome")
          // Module: prep_gtf
          PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, params.ref_gtf)
        } else {
          exit 1, "Pls supply --ref_gtf."
        }
      } else if (params.ref_fasta_ensembl) {
        // if ensembl name supplied:
        // Module: download ensembl genome
        DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // Module: prep_genome
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // Module: download ensembl gtf
        DOWNLOAD_FROM_ENSEMBL_GTF (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // Module: prep_gtf
        PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_ENSEMBL_GTF.out.gtf)
      } else if (params.ref_fasta_ucsc) {
        // Module: download ucsc genome
        DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // Module: download ucsc gtf
        DOWNLOAD_FROM_UCSC_GTF (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // Module: prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC_GTF.out.gtf)
        // Module: prep_gtf
        PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_UCSC_GTF.out.gtf)
      }
      // Module: prepare chromap index
      CHROMAP_INDEX (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name)
      // Module: run chromap atac
      CHROMAP_ATAC (MATCH_SAMPLE_NAME.out.sample_name.unique(), MATCH_SAMPLE_NAME.out.sample_files.collect(), CHROMAP_INDEX.out.index.first())
    }

    // GET_VALID_BARCODE
    //
    FILTER_CELL (sample_name_bam_fragment_valid_barcodes)
    // tuple val(sample_name), path(bam), path(fragment), path(filtered_barcode)


    // Emit PREP_GENOME output if PREP_GENOME is invoked.
    prep_genome_name    = Channel.empty()
    prep_genome_fasta   = Channel.empty()
    prep_gtf_genome     = Channel.empty()
    prep_gtf_file       = Channel.empty()

    try {
      prep_genome_name  = PREP_GENOME.out.genome_name
      // prep_genome       = Channel.value("run")
    } catch (Exception ex) {}
    try {
      prep_genome_fasta = PREP_GENOME.out.genome_fasta
    } catch (Exception ex) { }
    try {
      prep_gtf_genome   = PREP_GTF.out.genome_name
    } catch (Exception ex) {}
    try {
      prep_gtf_file = PREP_GTF.out.gtf
    } catch (Exception ex) { }

    res_files = Channel.empty()

  emit:
    res_files // out[0]: res folders for MultiQC report
    CHROMAP_ATAC.out.fragments // out[1]: for split bed
    CHROMAP_ATAC.out.ch_fragment // out[2]: fragment ch for ArchR
    "BAM_token1" // COMBINE_BAM.out.sample_name // out[3]: for split bam
    "BAM_token2" // COMBINE_BAM.out.bam // out[4]: for split bam
    prep_genome_name                      // out[5]: for downstream ArchR
    prep_genome_fasta                     // out[6]: for downstream ArchR
    prep_gtf_genome                       // out[7]: for downstream ArchR
    prep_gtf_file                         // out[8]: for downstream ArchR
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
