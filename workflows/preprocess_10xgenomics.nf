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
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['csv':'']] )
include { MATCH_SAMPLE_NAME } from '../modules/local/match_sample_name'
include { CELLRANGER_ATAC_COUNT } from '../modules/local/cellranger_atac_count'   addParams( options: modules['cellranger_atac_count'] )
include { DOWNLOAD_FROM_UCSC } from '../modules/local/download_from_ucsc'    addParams( options: modules['download_from_ucsc'] )
include { DOWNLOAD_FROM_ENSEMBL } from '../modules/local/download_from_ensembl'    addParams( options: modules['download_from_ensembl'] )
include { PREP_GENOME } from '../modules/local/prep_genome'
include { PREP_GTF } from '../modules/local/prep_gtf'
include { FILTER_CELL } from '../modules/local/filter_cell'
include { DOWNLOAD_FROM_UCSC_GTF } from '../modules/local/download_from_ucsc_gtf'    addParams( options: modules['download_from_ucsc_gtf'] )
include { DOWNLOAD_FROM_ENSEMBL_GTF } from '../modules/local/download_from_ensembl_gtf'    addParams( options: modules['download_from_ensembl_gtf'] )
include { CELLRANGER_INDEX } from '../modules/local/cellranger_index'             addParams( options: modules['cellranger_index'] )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow PREPROCESS_10XGENOMICS {
  take:
    reads
    sample_count

  main:
    // Examine if all required parameters supplied:
    if (!params.ref_cellranger_index && !(params.ref_fasta && params.ref_gtf) && !params.ref_fasta_ensembl && !params.ref_fasta_ucsc) {
      msg = "Must supply one from below:\n" + "Option1:\n  --ref_cellranger_index\n" + "Option2:\n  --ref_fasta\n  --ref_gtf\n" + "Option3:\n  --ref_fasta_ensembl\n" + "Option4:\n  --ref_fasta_ucsc\n"
      log.error msg
      exit 1, "EXIT!"
    }
    // Above is redundant to WorkflowMain::initialise()
    // module: staging sample_name in case of inconsistency in sample names
    MATCH_SAMPLE_NAME (reads, sample_count)

    if (params.ref_cellranger_index) {
      // if cellranger index folder provided:
      log.info "Parameter --ref_cellranger_index supplied, will use it as index folder."
      CELLRANGER_ATAC_COUNT (MATCH_SAMPLE_NAME.out.sample_name.unique(), MATCH_SAMPLE_NAME.out.sample_files.collect(), params.ref_cellranger_index)
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
      } else {
        exit 1, "PREPROCESS_10XGENOMICS: --ref_fasta_ucsc, or --ref_fasta_ensembl, or --ref_fasta/ref_gtf must be specified!"
      }
      // Module: prepare cellranger index
      CELLRANGER_INDEX (PREP_GENOME.out.genome_fasta, PREP_GTF.out.gtf, PREP_GENOME.out.genome_name)
      // Module: run cellranger-atac count
      CELLRANGER_ATAC_COUNT (MATCH_SAMPLE_NAME.out.sample_name.unique(), MATCH_SAMPLE_NAME.out.sample_files.collect(), CELLRANGER_INDEX.out.index_folder.first())
    }

    // Filter raw fragment and bam file based on filtered_peak_bc_matrix/barcodes.tsv
    FILTER_CELL (CELLRANGER_ATAC_COUNT.out.sample, "CR")

    // Emit PREP_GENOME output if PREP_GENOME is invoked.
    // prep_genome         = Channel.value("not_run")
    prep_genome_name    = Channel.empty()
    prep_genome_fasta   = Channel.empty()
    // prep_gtf            = Channel.value("not_run")
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
    res_files                             // out[0]: res folders for MultiQC report
    FILTER_CELL.out.filtered_fragment     // out[1]: for split bed
    FILTER_CELL.out.sample_name_filtered_fragment  // out[2]: fragment ch for ArchR
    FILTER_CELL.out.sample_name           // out[3]: for split bam
    FILTER_CELL.out.filtered_bam          // out[4]: for split bam
    prep_genome_name                      // out[5]: for downstream ArchR
    prep_genome_fasta                     // out[6]: for downstream ArchR
    prep_gtf_genome                       // out[7]: for downstream ArchR
    prep_gtf_file                         // out[8]: for downstream ArchR
}

// workflow.onComplete {
//     Completion.summary(workflow, params, log)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
