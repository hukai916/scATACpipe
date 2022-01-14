/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
// NfcoreSchema pops compilation error.

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
include { DOWNLOAD_FROM_UCSC } from '../modules/local/download_from_ucsc'    addParams( options: modules['download_from_ucsc'] )
include { DOWNLOAD_FROM_ENSEMBL } from '../modules/local/download_from_ensembl'    addParams( options: modules['download_from_ensembl'] )
include { BUILD_BSGENOME } from '../modules/local/build_bsgenome'
include { BUILD_TXDB } from '../modules/local/build_txdb'
include { PREP_GENOME } from '../modules/local/prep_genome'
include { PREP_GTF } from '../modules/local/prep_gtf'
include { BUILD_GENE_ANNOTATION } from '../modules/local/build_gene_annotation' addParams( options: modules['build_gene_annotation'] )
include { BUILD_GENOME_ANNOTATION } from '../modules/local/build_genome_annotation' addParams( options: modules['build_genome_annotation'] )
include { PREP_FRAGMENT } from '../modules/local/prep_fragment'
include { DOWNLOAD_FROM_UCSC_GTF } from '../modules/local/download_from_ucsc_gtf'    addParams( options: modules['download_from_ucsc_gtf'] )
include { DOWNLOAD_FROM_ENSEMBL_GTF } from '../modules/local/download_from_ensembl_gtf'    addParams( options: modules['download_from_ensembl_gtf'] )
include { AMULET_DETECT_DOUBLETS } from '../modules/local/amulet_detect_doublets'    addParams( options: modules['amulet_detect_doublets'] )
include { AMULET_MERGE_DOUBLETS } from '../modules/local/amulet_merge_doublets'
// For ArchR functions:
include { ARCHR_GET_ANNOTATION_BIOC } from '../modules/local/archr_get_annotation_bioc' addParams( options: modules['archr_get_annotation_bioc'] )
include { ARCHR_CREATE_ARROWFILES } from '../modules/local/archr_create_arrowfiles' addParams( options: modules['archr_create_arrowfiles'] )
include { ARCHR_CREATE_ARROWFILES_ANNOTATION } from '../modules/local/archr_create_arrowfiles_annotation' addParams( options: modules['archr_create_arrowfiles_annotation'] )
include { ARCHR_ADD_DOUBLETSCORES } from '../modules/local/archr_add_doubletscores' addParams( options: modules['archr_add_doubletscores'] )
include { ARCHR_ARCHRPROJECT } from '../modules/local/archr_archrproject' addParams( options: modules['archr_archrproject'] )
include { ARCHR_ARCHRPROJECT_ANNOTATION } from '../modules/local/archr_archrproject_annotation' addParams( options: modules['archr_archrproject_annotation'] )
include { ARCHR_FILTER_DOUBLETS } from '../modules/local/archr_filter_doublets' addParams( options: modules['archr_filter_doublets'] )
include { ARCHR_ARCHRPROJECT_QC } from '../modules/local/archr_archrproject_qc' addParams( options: modules['archr_archrproject_qc'] )
include { ARCHR_DIMENSION_REDUCTION } from '../modules/local/archr_dimension_reduction' addParams( options: modules['archr_dimension_reduction'] )
include { ARCHR_BATCH_CORRECTION } from '../modules/local/archr_batch_correction' addParams( options: modules['archr_batch_correction'] )
include { ARCHR_CLUSTERING } from '../modules/local/archr_clustering' addParams( options: modules['archr_clustering'] )
include { ARCHR_EMBEDDING } from '../modules/local/archr_embedding' addParams( options: modules['archr_embedding'] )
include { ARCHR_MARKER_GENE } from '../modules/local/archr_marker_gene' addParams( options: modules['archr_marker_gene'] )
include { ARCHR_SCRNASEQ_UNCONSTRAINED } from '../modules/local/archr_scrnaseq_unconstrained' addParams( options: modules['archr_scrnaseq_unconstrained'] )
include { ARCHR_SCRNASEQ_CONSTRAINED } from '../modules/local/archr_scrnaseq_constrained' addParams( options: modules['archr_scrnaseq_constrained'] )
include { ARCHR_PSEUDO_BULK_CLUSTERS } from '../modules/local/archr_pseudo_bulk_clusters' addParams( options: modules['archr_pseudo_bulk_clusters'] )
include { ARCHR_PSEUDO_BULK_CLUSTERS2 } from '../modules/local/archr_pseudo_bulk_clusters2' addParams( options: modules['archr_pseudo_bulk_clusters2'] )
include { ARCHR_CALL_PEAKS_CLUSTERS } from '../modules/local/archr_call_peaks_clusters' addParams( options: modules['archr_call_peaks_clusters'] )
include { ARCHR_CALL_PEAKS_CLUSTERS2 } from '../modules/local/archr_call_peaks_clusters2' addParams( options: modules['archr_call_peaks_clusters2'] )
include { ARCHR_GET_MARKER_PEAKS_CLUSTERS } from '../modules/local/archr_get_marker_peaks_clusters' addParams( options: modules['archr_get_marker_peaks_clusters'] )
include { ARCHR_GET_MARKER_PEAKS_CLUSTERS2 } from '../modules/local/archr_get_marker_peaks_clusters2' addParams( options: modules['archr_get_marker_peaks_clusters2'] )
include { ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS } from '../modules/local/archr_marker_peaks_in_tracks_clusters' addParams( options: modules['archr_marker_peaks_in_tracks_clusters'] )
include { ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2 } from '../modules/local/archr_marker_peaks_in_tracks_clusters2' addParams( options: modules['archr_marker_peaks_in_tracks_clusters2'] )
include { ARCHR_PAIRWISE_TEST_CLUSTERS } from '../modules/local/archr_pairwise_test_clusters' addParams( options: modules['archr_pairwise_test_clusters'] )
include { ARCHR_PAIRWISE_TEST_CLUSTERS2 } from '../modules/local/archr_pairwise_test_clusters2' addParams( options: modules['archr_pairwise_test_clusters2'] )
include { ARCHR_MOTIF_ENRICHMENT_CLUSTERS } from '../modules/local/archr_motif_enrichment_clusters' addParams( options: modules['archr_motif_enrichment_clusters'] )
include { ARCHR_MOTIF_ENRICHMENT_CLUSTERS2 } from '../modules/local/archr_motif_enrichment_clusters2' addParams( options: modules['archr_motif_enrichment_clusters2'] )
include { ARCHR_MOTIF_DEVIATIONS_CLUSTERS } from '../modules/local/archr_motif_deviations_clusters' addParams( options: modules['archr_motif_deviations_clusters'] )
include { ARCHR_MOTIF_DEVIATIONS_CLUSTERS2 } from '../modules/local/archr_motif_deviations_clusters2' addParams( options: modules['archr_motif_deviations_clusters2'] )
include { ARCHR_FOOTPRINTING_CLUSTERS } from '../modules/local/archr_footprinting_clusters' addParams( options: modules['archr_footprinting_clusters'] )
include { ARCHR_FOOTPRINTING_CLUSTERS2 } from '../modules/local/archr_footprinting_clusters2' addParams( options: modules['archr_footprinting_clusters2'] )
include { ARCHR_COACCESSIBILITY_CLUSTERS } from '../modules/local/archr_coaccessibility_clusters' addParams( options: modules['archr_coaccessibility_clusters'] )
include { ARCHR_COACCESSIBILITY_CLUSTERS2 } from '../modules/local/archr_coaccessibility_clusters2' addParams( options: modules['archr_coaccessibility_clusters2'] )
include { ARCHR_PEAK2GENELINKAGE_CLUSTERS2 } from '../modules/local/archr_peak2genelinkage_clusters2' addParams( options: modules['archr_peak2genelinkage_clusters2'] )
include { ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS } from '../modules/local/archr_get_positive_tf_regulator_clusters' addParams( options: modules['archr_get_positive_tf_regulator_clusters'] )
include { ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2 } from '../modules/local/archr_get_positive_tf_regulator_clusters2' addParams( options: modules['archr_get_positive_tf_regulator_clusters2'] )
include { ARCHR_TRAJECTORY_CLUSTERS2 } from '../modules/local/archr_trajectory_clusters2' addParams( options: modules['archr_trajectory_clusters2'] )
include { ARCHR_GET_CLUSTERING_TSV } from '../modules/local/archr_get_clustering_tsv' addParams( options: modules['archr_get_clustering_tsv'] )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow DOWNSTREAM_ARCHR {
  take:
    fragments // tuple val(sample_name), path(fragment)
    with_preprocess // string value: "preprocess_null", "preprocess_default", or "preprocess_10xgenomics"
    prep_genome // string value: "run" or "not_run"
    prep_genome_name // PREP_GENOME.out.genome_name if prep_genome == 'run' else Channel.empty()
    prep_genome_fasta // PREP_GENOME.out.genome_fasta if prep_genome == 'run' else Channel.empty()
    prep_gtf   // string value: "run" or "not_run"
    prep_gtf_genome // PREP_GTF.out.genome_name if prep_gtf == 'run' else Channel.empty()
    prep_gtf_file // PREP_GTF.out.gtf if prep_gtf == 'run' else Channel.empty()
    // Won't take tuple as input

  main:
    // Examine if all required parameters supplied:
    def archr_input_type = "" // "naive", "genome_gtf", "bsgenome_txdb_org"
    def archr_input_list = []
    // for "genome_gtf": [genome_name, genome_fasta, gtf]
    // for "bsgenome_txdb_org"" [bsgenome, txdb, org]
    if (params.archr_bsgenome && params.archr_org && params.archr_txdb) {
      archr_input_type = "bsgenome_txdb_org"
    } else if (params.archr_genome_fasta && params.ref_gtf) {
      archr_input_type = "genome_gtf"
    } else if (with_preprocess == "preprocess_null") {
      if (!params.archr_genome) {
        msg = "ArchR genome required but not supplied!\nOption1:\n  --archr_genome [a genome name]\nOption2:\n  --archr_genome_fasta [path to genome fasta]\n  --ref_gtf [path to gtf file]\n  --archr_blacklist [optional, path to blacklist file]\nOption3:\n  --archr_bsgenome [Bioconductor BSgenome name]\n  --archr_txdb [Bioconductor TxDb name]\n  --archr_org [Bioconductor OrgDb name]\n  --archr_blacklist [optional, path to blacklist file]\nPlease supply the above params to continue.\n"
        log.error msg
        exit 1, "EXIT!"
      }
    } else if ((with_preprocess == "preprocess_default") || (with_preprocess == "preprocess_10xgenomics") || (with_preprocess == "preprocess_chromap")) {
      msg = "ArchR genome required but not supplied!\nOption1:\n  --ref_fasta_ucsc [a genome name]\nOption2:\n  --ref_fasta_ensembl [a genome name]\nOption3:\n  --ref_fasta [path to genome fasta]\n  --ref_gtf [path to gtf file]\nOption4:\n  --archr_genome_fasta [path to genome fasta]\n  --ref_gtf [path to gtf file]\n  --archr_blacklist [optional, path to blacklist file]\nOption5:\n  --archr_bsgenome [path to BSgenome obj]\n  --archr_txdb [path to TxDb obj]\n  --archr_org [path to OrgDb obj]\n  --archr_blacklist [optional, path to blacklist file]\nPlease supply the above params to continue.\n"
      if (params.ref_bwa_index || params.ref_cellranger_index) {
        if (!params.ref_fasta_ensembl && !params.ref_fasta_ucsc) {
          log.error msg
          exit 1, "EXIT!"
        }
      } else if (params.ref_fasta) {
        if (!params.ref_gtf) {
          log.error msg
          exit 1, "EXIT!"
        }
      } else if (!params.ref_fasta_ensembl && !params.ref_fasta_ucsc) {
        log.error msg
        exit 1, "EXIT!"
      }
    }
    // Above is redundant to WorkflowMain::initialise()

    if (params.archr_bsgenome && params.archr_org && params.archr_txdb) {
      log.info "INFO: --archr_bsgenome, --archr_txdb, and --archr_org supplied."
      archr_input_type = "bsgenome_txdb_org"
      archr_input_list = [params.archr_bsgenome, params.archr_txdb, params.archr_org]
    } else if (params.archr_genome_fasta && params.ref_gtf) {
      log.info "INFO: --archr_genome_fasta and --ref_gtf supplied."
      PREP_GENOME(Channel.fromPath(params.test_fasta), "custom_genome")
      PREP_GTF(PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, params.ref_gtf)
      archr_input_type = "genome_gtf"
      archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
    } else {
      if (with_preprocess == "preprocess_null") {
        if (params.archr_genome) {
          log.info "INFO: --archr_genome supplied."

          if (["hg38", "hg19", "mm10", "mm9"].contains(params.archr_genome)) {
            log.info "INFO: natively supported ArchR genome supplied."

            archr_input_type = "naive"
            archr_input_list = [params.archr_genome, "NA", "NA"]
          } else if (genome_ensembl_list.contains(params.archr_genome)) {
              DOWNLOAD_FROM_ENSEMBL (params.archr_genome, Channel.fromPath('assets/genome_ensembl.json'))
              DOWNLOAD_FROM_ENSEMBL_GTF (params.archr_genome, Channel.fromPath('assets/genome_ensembl.json'))
              PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
              PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_ENSEMBL_GTF.out.gtf)
              archr_input_type = "genome_gtf"
              archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
          } else if (genome_ucsc_list.contains(params.archr_genome)) {
              DOWNLOAD_FROM_UCSC (params.archr_genome, Channel.fromPath('assets/genome_ucsc.json'))
              DOWNLOAD_FROM_UCSC_GTF(params.archr_genome, Channel.fromPath('assets/genome_ucsc.json'))
              PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_name)
              PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_UCSC_GTF.out.gtf)
              archr_input_type = "genome_gtf"
              archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
          } else {
              exit 1, "Pls use the --support_genome to show a list of supported genomes!"
          }
        } else {
          exit 1, "ArchR genome required but not supplied!\nOption1:\n  --archr_genome [a genome name]\nOption2:\n  --archr_genome_fasta [path to genome fasta]\n  --ref_gtf [path to gtf file]\n  --archr_blacklist [optional, path to blacklist file]\nOption3:\n  --archr_bsgenome [Bioconductor BSgenome name]\n  --archr_txdb [Bioconductor TxDb name]\n  --archr_org [Bioconductor OrgDb name]\n  --archr_blacklist [optional, path to blacklist file]\nPlease supply the above params to continue.\n"
        }
      } else if (with_preprocess == "preprocess_default" || with_preprocess == "preprocess_10xgenomics" || with_preprocess == "preprocess_chromap") {
        exit_msg = "ArchR genome required but not supplied!\nOption1:\n  --ref_fasta_ucsc [a genome name]\nOption2:\n  --ref_fasta_ensembl [a genome name]\nOption3:\n  --ref_fasta [path to genome fasta]\n  --ref_gtf [path to gtf file]\nOption4:\n  --archr_genome_fasta [path to genome fasta]\n  --ref_gtf [path to gtf file]\n  --archr_blacklist [optional, path to blacklist file]\nOption5:\n  --archr_bsgenome [path to BSgenome obj]\n  --archr_txdb [path to TxDb obj]\n  --archr_org [path to OrgDb obj]\n  --archr_blacklist [optional, path to blacklist file]\nPlease supply the above params to continue.\n"
        if (params.ref_bwa_index || params.ref_cellranger_index) {
          // Need to download genome and gtf:
          if (params.ref_fasta_ensembl) {
            DOWNLOAD_FROM_ENSEMBL(params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
            DOWNLOAD_FROM_ENSEMBL_GTF(params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
            PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
            PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_ENSEMBL_GTF.out.gtf)
            archr_input_type = "genome_gtf"
            archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
          } else if (params.ref_fasta_ucsc) {
            DOWNLOAD_FROM_UCSC(params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
            DOWNLOAD_FROM_UCSC_GTF(params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
            PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_name)
            PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_UCSC_GTF.out.gtf)
            archr_input_type = "genome_gtf"
            archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
          } else {
            exit 1, exit_msg
          }
        } else if (params.ref_fasta) {
          // PREP_GENOME should has been performed
          // If PREPROCESS_DEFAULT: PREP_GTF has not been performed
          // If PREPROCESS_10XGENOMICS: PREP_GTF has been performed
          if (params.ref_gtf) {
            tem_genome_name  = Channel.empty()
            tem_genome_fasta = Channel.empty()
            tem_gtf_file     = Channel.empty()

            if (!(prep_genome == "run")) {
              log.error 'Something went wrong!'
              exit 1, "Exit!"
            } else {
              tem_genome_name = prep_genome_name
              tem_genome_fasta = prep_genome_fasta
            }

            if (prep_gtf == "run") {
              tem_gtf_file = prep_gtf_file
            } else {
              PREP_GTF (prep_genome_fasta, prep_genome_name, params.ref_gtf)
              tem_gtf_file = PREP_GTF.out.gtf
            }

            archr_input_type = "genome_gtf"
            archr_input_list = [tem_genome_name.collect(), tem_genome_fasta.collect(), tem_gtf_file.collect()]
          } else {
            exit 1, "Pls also supply --ref_gtf."
          }
        } else if (params.ref_fasta_ensembl) {
          if (with_preprocess == "preprocess_default") {
            // If PREPROCESS_DEFAULT: PREP_GENOME has been performed but not PREP_GTF
            if (!(prep_genome == "run")) {
              log.error "Something must be wrong!"
              exit 1, "EXIT!"
            }
            DOWNLOAD_FROM_ENSEMBL_GTF(params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
            PREP_GTF (prep_genome_name, prep_genome_fasta, DOWNLOAD_FROM_ENSEMBL_GTF.out.gtf)
            archr_input_type = "genome_gtf"
            // archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
            archr_input_list = [prep_genome_name.collect(), prep_genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
          } else if (with_preprocess == "preprocess_10xgenomics" || with_preprocess == "preprocess_chromap") {
            // If PREPPROCESS_10XGENOMICS: both PREP_GENOME and PREP_GTF should been performed
            if (!(prep_genome == "run") || !(prep_gtf == "run")) {
              log.error "Something must be wrong(2)!"
              exit 1, "EXIT!"
            }
            archr_input_type = "genome_gtf"
            // archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
            archr_input_list = [prep_genome_name.collect(), prep_genome_fasta.collect(), prep_gtf_file.collect()]
          }
        } else if (params.ref_fasta_ucsc) {
          if (with_preprocess == "preprocess_default") {
            // If PREPROCESS_DEFAULT: PREP_GENOME has been performed but not PREP_GTF
            if (!(prep_genome == "run")) {
              log.error "Something must be wrong!"
              exit 1, "EXIT!"
            }
            // Check if natively supported for ucsc genomes:
            if (["hg38", "hg19", "mm10", "mm9"].contains(params.ref_fasta_ucsc)) {
              log.info "INFO: natively supported ArchR genome supplied."
              archr_input_type = "naive"
              archr_input_list = [params.ref_fasta_ucsc, "NA", "NA"]
            } else {
              DOWNLOAD_FROM_UCSC_GTF(params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
              // PREP_GTF (PREP_GENOME.out.genome_fasta, PREP_GENOME.out.genome_name, DOWNLOAD_FROM_UCSC_GTF.out.gtf)
              PREP_GTF (prep_genome_fasta, prep_genome_name, DOWNLOAD_FROM_UCSC_GTF.out.gtf)
              archr_input_type = "genome_gtf"
              // archr_input_list = [PREP_GENOME.out.genome_name.collect(), PREP_GENOME.out.genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
              archr_input_list = [prep_genome_name.collect(), prep_genome_fasta.collect(), PREP_GTF.out.gtf.collect()]
            }
          } else if (with_preprocess == "preprocess_10xgenomics" || with_preprocess == "preprocess_chromap") {
            // If PREPPROCESS_10XGENOMICS: both PREP_GENOME and PREP_GTF should been performed
            if (!(prep_genome == "run") || !(prep_gtf == "run")) {
              log.error "Something must be wrong (2)!"
              exit 1, "EXIT!"
            }
            // Check if natively supported for ucsc genomes:
            if (["hg38", "hg19", "mm10", "mm9"].contains(params.ref_fasta_ucsc)) {
              log.info "INFO: natively supported ArchR genome supplied."
              archr_input_type = "naive"
              archr_input_list = [params.ref_fasta_ucsc, "NA", "NA"]
            } else {
              archr_input_type = "genome_gtf"
              archr_input_list = [prep_genome_name.collect(), prep_genome_fasta.collect(), prep_gtf_file.collect()]
            }
          }
        } else {
          exit 1, exit_msg
        }
      }
    }

    log.info "archr_input_type: " + archr_input_type
    // Depending on ArchR input type, prepare ArchR annotation files accordingly:
    if (archr_input_type == "naive") {
      // Run ArchR normally:
      log.info "Naively supported ArchR genome: " + archr_input_list[0] + " will be used."

      ARCHR_CREATE_ARROWFILES(fragments, archr_input_list[0], params.archr_thread)
      // Module: add DoubletScores
      ARCHR_ADD_DOUBLETSCORES(ARCHR_CREATE_ARROWFILES.out.sample_name, ARCHR_CREATE_ARROWFILES.out.arrowfile)
      // ch_samplename_list = ARCHR_ADD_DOUBLETSCORES.out.sample_name.toSortedList()
      ch_arrowfile_list = ARCHR_ADD_DOUBLETSCORES.out.arrowfile.toSortedList( { a, b -> a.getName() <=> b.getName() })
      ARCHR_ARCHRPROJECT(ch_arrowfile_list, archr_input_list[0], params.archr_thread)
      ARCHR_ARCHRPROJECT_QC(ARCHR_ARCHRPROJECT.out.archr_project)
    } else if (archr_input_type == "bsgenome_txdb_org") {
      // Note that for this option, all supplied package names must be available from Bioconductor per .requirePackage() requirement.
      // Run ArchR with ANNOTATION option
      log.info "INFO: ArchR will build gene/genomeAnnotation files with custom TxDb, Org, and BSgenome files supplied by user."

      ARCHR_GET_ANNOTATION_BIOC(params.archr_txdb, params.archr_org, params.archr_bsgenome)
      ARCHR_CREATE_ARROWFILES_ANNOTATION(fragments, ARCHR_GET_ANNOTATION_BIOC.out.geneAnnotation.collect(), ARCHR_GET_ANNOTATION_BIOC.out.genomeAnnotation.collect(), ARCHR_GET_ANNOTATION_BIOC.out.user_rlib.collect(), params.archr_thread)
      // Module: add DoubletScores
      ARCHR_ADD_DOUBLETSCORES(ARCHR_CREATE_ARROWFILES_ANNOTATION.out.sample_name, ARCHR_CREATE_ARROWFILES_ANNOTATION.out.arrowfile)
      // ch_samplename_list = ARCHR_ADD_DOUBLETSCORES.out.sample_name.toSortedList()
      ch_arrowfile_list = ARCHR_ADD_DOUBLETSCORES.out.arrowfile.toSortedList( { a, b -> a.getName() <=> b.getName() })

      ARCHR_ARCHRPROJECT_ANNOTATION(ch_arrowfile_list, ARCHR_GET_ANNOTATION_BIOC.out.geneAnnotation, ARCHR_GET_ANNOTATION_BIOC.out.genomeAnnotation, ARCHR_GET_ANNOTATION_BIOC.out.user_rlib)
      ARCHR_ARCHRPROJECT_QC(ARCHR_ARCHRPROJECT_ANNOTATION.out.archr_project)
    } else if (archr_input_type == "genome_gtf") {
      if (params.species_latin_name) {
        // Build BSgenome:
        if (prep_genome == "run") {
          BUILD_BSGENOME(prep_genome_fasta)
        } else if (prep_genome == "not_run") {
          BUILD_BSGENOME(PREP_GENOME.out.genome_fasta)
        }

        // Build ArchR gene annotation file:
        if (prep_gtf == "run") {
          BUILD_TXDB (BUILD_BSGENOME.out.bsgenome, prep_gtf_file)
          BUILD_GENE_ANNOTATION(BUILD_TXDB.out.txdb, prep_gtf_file, params.species_latin_name)
        } else if (prep_gtf == "not_run") {
          BUILD_TXDB (BUILD_BSGENOME.out.bsgenome, PREP_GTF.out.gtf)
          BUILD_GENE_ANNOTATION(BUILD_TXDB.out.txdb, PREP_GTF.out.gtf, params.species_latin_name)
        }

        // Build ArchR genome annotation file:
        if (params.archr_blacklist) {
          BUILD_GENOME_ANNOTATION(BUILD_BSGENOME.out.bsgenome, BUILD_GENE_ANNOTATION.out.gene_annotation, Channel.fromPath(params.archr_blacklist))
        } else {
          BUILD_GENOME_ANNOTATION(BUILD_BSGENOME.out.bsgenome, BUILD_GENE_ANNOTATION.out.gene_annotation, "$projectDir/assets/file_token.txt")
        }

        // Match fragment file against gtf file:
        PREP_FRAGMENT(fragments, archr_input_list[2])
        ARCHR_CREATE_ARROWFILES_ANNOTATION(PREP_FRAGMENT.out.fragments, BUILD_GENE_ANNOTATION.out.gene_annotation.collect(), BUILD_GENOME_ANNOTATION.out.genome_annotation.collect(), BUILD_BSGENOME.out.user_rlib.collect(), params.archr_thread)
        // Module: add DoubletScores
        ARCHR_ADD_DOUBLETSCORES(ARCHR_CREATE_ARROWFILES_ANNOTATION.out.sample_name, ARCHR_CREATE_ARROWFILES_ANNOTATION.out.arrowfile)
        // ch_samplename_list = ARCHR_ADD_DOUBLETSCORES.out.sample_name.toSortedList()
        ch_arrowfile_list = ARCHR_ADD_DOUBLETSCORES.out.arrowfile.toSortedList( { a, b -> a.getName() <=> b.getName() })

        ARCHR_ARCHRPROJECT_ANNOTATION(ch_arrowfile_list, BUILD_GENE_ANNOTATION.out.gene_annotation, BUILD_GENOME_ANNOTATION.out.genome_annotation, BUILD_BSGENOME.out.user_rlib)
        ARCHR_ARCHRPROJECT_QC(ARCHR_ARCHRPROJECT_ANNOTATION.out.archr_project)
      } else {
        exit 1, "Pls also supply --species_latin_name."
      }
    }

    // Module: AMULET doublet filtering
    if (params.doublet_removal_algorithm == "amulet") {
      if (archr_input_type == "genome_gtf") {
        // Use prep_fragment.out.fragments
        AMULET_DETECT_DOUBLETS(PREP_FRAGMENT.out.fragments, Channel.fromPath(params.amulet_rmsk_bed).first(), Channel.fromPath(params.amulet_autosomes).first())
      } else {
        // Use fragments
        AMULET_DETECT_DOUBLETS(fragments, Channel.fromPath(params.amulet_rmsk_bed).first(), Channel.fromPath(params.amulet_autosomes).first())
      }
      // Module: generate Doublet cell list input to ArchR
      AMULET_MERGE_DOUBLETS(AMULET_DETECT_DOUBLETS.out.cells_filter.collect())
    }

    // Module: filter doublets depending on user option.
    if (!params.archr_filter_doublets_ratio) {
      // Module: dimension reduction
      ARCHR_DIMENSION_REDUCTION(ARCHR_ARCHRPROJECT_QC.out.archr_project)
    } else if (params.doublet_removal_algorithm == "archr") { // for test only
      // Module: filtering doublets
      ARCHR_FILTER_DOUBLETS(ARCHR_ARCHRPROJECT_QC.out.archr_project, params.archr_filter_doublets_ratio)
      // Module: dimension reduction
      ARCHR_DIMENSION_REDUCTION(ARCHR_FILTER_DOUBLETS.out.archr_project)
    } else if (params.doublet_removal_algorithm == "amulet") {
      // Module: filtering doublets
      AMULET_FILTER_DOUBLETS(ARCHR_ARCHRPROJECT_QC.out.archr_project, AMULET_DETECT_DOUBLETS.out.cells_filter)
      // Module: dimension reduction
      ARCHR_DIMENSION_REDUCTION(AMULET_FILTER_DOUBLETS.out.archr_project)
    }

    // Module: batch correction with harmony
    ARCHR_BATCH_CORRECTION(ARCHR_DIMENSION_REDUCTION.out.archr_project)

    // Module: clustering with Seurat's FindClusters() function
    ARCHR_CLUSTERING(ARCHR_BATCH_CORRECTION.out.archr_project)

    // Module: single-cell embeddings
    ARCHR_EMBEDDING(ARCHR_CLUSTERING.out.archr_project)

    // Module: find marker gene
    ARCHR_MARKER_GENE(ARCHR_EMBEDDING.out.archr_project)

    // Module: integrate with matching scRNAseq data
    if (!(params.archr_scrnaseq)) {
      params.groupby_cluster = "Clusters"
      log.info "NOTICE: --archr_scrnaseq: not supplied, skip integrative analysis with scRNA-seq!"
      // ARCHR_PSEUDO_BULK(ARCHR_MARKER_GENE.out.archr_project, params.groupby_cluster)
      ARCHR_PSEUDO_BULK_CLUSTERS(ARCHR_MARKER_GENE.out.archr_project)
      // For each Arrorproject, you can have only one set of peak set unless you copy arrow files and create another arrowproject. That is why we implemented ARCHR_PSEUDO_BULK_CLUSTERS and ARCHR_PSEUDO_BULK_CLUSTERS2
    } else {
        params.groupby_cluster = "Clusters2"
        log.info "NOTICE: --archr_scrnaseq: supplied, will perform integrative analysis with scRNA-seq!"
        ARCHR_PSEUDO_BULK_CLUSTERS(ARCHR_MARKER_GENE.out.archr_project)
        ARCHR_SCRNASEQ_UNCONSTRAINED(ARCHR_MARKER_GENE.out.archr_project, params.archr_scrnaseq)
        // log.info "INFO: use the following cluster names to define --archr_scrnaseq_grouplist."
        ARCHR_SCRNASEQ_UNCONSTRAINED.out.cell_type_scrna
          .splitText()
          .subscribe onNext: { String str -> println "Cluster name from scRNAseq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --archr_scrnaseq_grouplist ***\n"}

        if ((!params.archr_scrnaseq_grouplist)) {
          log.info "NOTICE: --archr_scrnaseq_grouplist: not supplied, skip constrained integration!"
          // ARCHR_PSEUDO_BULK(ARCHR_SCRNASEQ_UNCONSTRAINED.out.archr_project, params.groupby_cluster)
          ARCHR_PSEUDO_BULK_CLUSTERS2(ARCHR_SCRNASEQ_UNCONSTRAINED.out.archr_project)
        } else {
            log.info "NOTICE: --archr_scrnaseq_grouplist: supplied, will perform constrained integration!"
            ARCHR_SCRNASEQ_CONSTRAINED(ARCHR_SCRNASEQ_UNCONSTRAINED.out.archr_project, params.archr_scrnaseq, params.archr_scrnaseq_grouplist)
            // ARCHR_PSEUDO_BULK(ARCHR_SCRNASEQ_CONSTRAINED.out.archr_project, params.groupby_cluster)
            ARCHR_PSEUDO_BULK_CLUSTERS2(ARCHR_SCRNASEQ_CONSTRAINED.out.archr_project)
        }
    }

    // Module: call peaks
    if (params.groupby_cluster == "Clusters") {
      ARCHR_CALL_PEAKS_CLUSTERS(ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_project)
    } else if (params.groupby_cluster == "Clusters2") {
      ARCHR_CALL_PEAKS_CLUSTERS(ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_project)
      ARCHR_CALL_PEAKS_CLUSTERS2(ARCHR_PSEUDO_BULK_CLUSTERS2.out.archr_project)
    }

    // Module: identify marker peaks and perform MA/Volcano plots
    if (params.groupby_cluster == "Clusters") {
      ARCHR_GET_MARKER_PEAKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project)
      ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.group_names
        .splitText()
        .subscribe onNext: { String str -> println "Group name from scATAC-seq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --pairwise_test_clusters_1/2 and --marker_peak_clusters***\n"}
    } else if (params.groupby_cluster == "Clusters2") {
      ARCHR_GET_MARKER_PEAKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project)
      ARCHR_GET_MARKER_PEAKS_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project)
      ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.group_names
        .splitText()
        .subscribe onNext: { String str -> println "Group name from scATAC-seq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --pairwise_test_clusters_1/2 ***\n"}
      ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.group_names
        .splitText()
        .subscribe onNext: { String str -> println "Group name from scRNA-seq: ${str}".trim() }, onComplete: { print "\n*** use above names to define --pairwise_test_clusters2_1/2 and --marker_peak_clusters2***\n"}
    }

    // Module: plot peaks in browser tracks
    if (params.groupby_cluster == "Clusters") {
      if (!(params.marker_peak_geneSymbol && params.marker_peak_clusters)) {
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: not supplied, skip marker peak plotting on browser tracks!"
      } else {
        // Perform plotting
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: supplied, will perform marker peak plotting on browser tracks!"
        ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.marker_peak_geneSymbol, params.marker_peak_clusters)
      }
    } else if (params.groupby_cluster == "Clusters2") {
      if (!(params.marker_peak_geneSymbol && params.marker_peak_clusters)) {
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: not supplied, skip marker peak plotting on browser tracks!"
      } else {
          // Perform plotting
          log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters: supplied, will perform marker peak plotting on browser tracks!"
          ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.marker_peak_geneSymbol, params.marker_peak_clusters)
      }

      if (!(params.marker_peak_geneSymbol && params.marker_peak_clusters2)) {
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters2: not supplied, skip marker peak plotting on browser tracks!"
      } else {
        // Perform plotting
        log.info "NOTICE: --marker_peak_geneSymbol and --marker_peak_clusters2: supplied, will perform marker peak plotting on browser tracks!"
        ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project, ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.archr_marker_peaks, params.marker_peak_geneSymbol, params.marker_peak_clusters2)
      }
    }

    // Module: perform pairwise test
    if (params.groupby_cluster == "Clusters") {
      if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
        log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip pairwise plotting!"
      } else {
        // Perform plotting
        log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform pairwise plotting!"
        ARCHR_PAIRWISE_TEST_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2)
      }
    } else if (params.groupby_cluster == "Clusters2") {
        if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
          log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip pairwise plotting!"
      } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform pairwise plotting!"
          ARCHR_PAIRWISE_TEST_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2)
      }

      if (!(params.pairwise_test_clusters2_1 && params.pairwise_test_clusters2_2)) {
        log.info "NOTICE: --pairwise_test_clusters2_1/2: not supplied, skip pairwise plotting!"
      } else {
        // Perform plotting
        log.info "NOTICE: --pairwise_test_clusters2_1/2: supplied, perform pairwise plotting!"
        ARCHR_PAIRWISE_TEST_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project, params.pairwise_test_clusters2_1, params.pairwise_test_clusters2_2)
      }
    }

    // Module: motif enrichment: note that ARCHR_MOTIF_ENRICHMENT_CLUSTERS and ARCHR_MOTIF_ENRICHMENT_CLUSTERS2 are exactly the same except for the outdir name.
    if (params.groupby_cluster == "Clusters") {
      if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
        log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip motif enrichment!"
      } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform motif enrichment!"
          ARCHR_MOTIF_ENRICHMENT_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_PAIRWISE_TEST_CLUSTERS.out.archr_marker_test, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2, params.custom_peaks)
      }
    } else if (params.groupby_cluster == "Clusters2") {
        if (!(params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2)) {
          log.info "NOTICE: --pairwise_test_clusters_1/2: not supplied, skip motif enrichment!"
        } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters_1/2: supplied, perform motif enrichment!"
          ARCHR_MOTIF_ENRICHMENT_CLUSTERS(ARCHR_CALL_PEAKS_CLUSTERS.out.archr_project, ARCHR_PAIRWISE_TEST_CLUSTERS.out.archr_marker_test, ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.archr_marker_peaks, params.pairwise_test_clusters_1, params.pairwise_test_clusters_2, params.custom_peaks)
        }

      if (!(params.pairwise_test_clusters2_1 && params.pairwise_test_clusters2_2)) {
        log.info "NOTICE: --pairwise_test_clusters2_1/2: not supplied, skip motif enrichment!"
      } else {
          // Perform plotting
          log.info "NOTICE: --pairwise_test_clusters2_1/2: supplied, perform motif enrichment!"
          ARCHR_MOTIF_ENRICHMENT_CLUSTERS2(ARCHR_CALL_PEAKS_CLUSTERS2.out.archr_project, ARCHR_PAIRWISE_TEST_CLUSTERS2.out.archr_marker_test, ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.archr_marker_peaks, params.pairwise_test_clusters2_1, params.pairwise_test_clusters2_2, params.custom_peaks)
      }
    }

    if ((params.pairwise_test_clusters_1 && params.pairwise_test_clusters_2) || (params.pairwise_test_clusters2_1 && params.pairwise_test_clusters2_2)) {
      // Module: motif deviation,require motif enrichment result
      if (params.groupby_cluster == "Clusters") {
        ARCHR_MOTIF_DEVIATIONS_CLUSTERS(ARCHR_MOTIF_ENRICHMENT_CLUSTERS.out.archr_project, params.custom_peaks)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_MOTIF_DEVIATIONS_CLUSTERS(ARCHR_MOTIF_ENRICHMENT_CLUSTERS.out.archr_project, params.custom_peaks)
          ARCHR_MOTIF_DEVIATIONS_CLUSTERS2(ARCHR_MOTIF_ENRICHMENT_CLUSTERS2.out.archr_project, params.custom_peaks)
      }

      // Module: footprinting
      if (params.groupby_cluster == "Clusters") {
        ARCHR_FOOTPRINTING_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project, ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_dir)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_FOOTPRINTING_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project, ARCHR_PSEUDO_BULK_CLUSTERS.out.archr_dir)
          ARCHR_FOOTPRINTING_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project, ARCHR_PSEUDO_BULK_CLUSTERS2.out.archr_dir)
      }

      // Module: integrative analysis
      // Below are for integrative analysis: co-accessibility; peak2genelinkage; positive TF regulators.
      // Module: co-accessibility (for both clusters and clusters2)
      if (params.groupby_cluster == "Clusters") {
        ARCHR_COACCESSIBILITY_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_COACCESSIBILITY_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
          ARCHR_COACCESSIBILITY_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project)
      }

      // Module: peak2genelinkage: for clusters2 only
      if (params.groupby_cluster == "Clusters2") {
        ARCHR_PEAK2GENELINKAGE_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project)
      }

      // Module: identify "positive" TF-regulators
      if (params.groupby_cluster == "Clusters") {
        ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
      } else if (params.groupby_cluster == "Clusters2") {
          ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.archr_project)
          ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project)
      }

      // Module: trajectory analysis: for Clusters2 only
      // TODO: Module: trajectory analysis: for Clusters using Gene Score Matrix
      if (params.groupby_cluster == "Clusters2") {
        if (!params.trajectory_groups) {
          log.info "Parameter --trajectory_groups not supplied, checking trajectory analysis!"
        } else {
            log.info "Parameter --trajectory_groups supplied, will perform trajectory analysis!"
            ARCHR_TRAJECTORY_CLUSTERS2(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.archr_project, params.trajectory_groups)
        }
      } else {
        log.info "Parameter --scrnaseq not supplied, skip trajectory analysis!"
      }
    }

    // Module: prepare clustering tsv file for spliting using sinto fragment
    if (params.groupby_cluster == "Clusters") {
      // ARCHR_CLUSTERING.out.archr_project.view()
      if (archr_input_type == "genome_gtf") {
        ARCHR_GET_CLUSTERING_TSV(ARCHR_CLUSTERING.out.archr_project.collect(), PREP_FRAGMENT.out.fragments, "Clusters")
      } else {
        ARCHR_GET_CLUSTERING_TSV(ARCHR_CLUSTERING.out.archr_project.collect(), fragments, "Clusters")
      }
    } else if (params.groupby_cluster == "Clusters2") {
      if (archr_input_type == "genome_gtf") {
        ARCHR_GET_CLUSTERING_TSV(ARCHR_PSEUDO_BULK_CLUSTERS2.out.archr_project.collect(), PREP_FRAGMENT.out.fragments, "Clusters2")
      } else {
        ARCHR_GET_CLUSTERING_TSV(ARCHR_PSEUDO_BULK_CLUSTERS2.out.archr_project.collect(), fragments, "Clusters2")
      }

    }

    // Collect all output results for MultiQC report:
    res_files = Channel.empty()

    // ARCHR_CREATE_ARROWFILES module:
    try {
      res_files = res_files.mix(ARCHR_CREATE_ARROWFILES.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_ADD_DOUBLETSCORES
    try {
      res_files = res_files.mix(ARCHR_ADD_DOUBLETSCORES.out.report.collect().ifEmpty([]))
      res_files = res_files.mix(ARCHR_ADD_DOUBLETSCORES.out.summary.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_ARCHRPROJECT_QC:
    try {
      res_files = res_files.mix(ARCHR_ARCHRPROJECT_QC.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_FILTER_DOUBLETS:
    try {
      res_files = res_files.mix(ARCHR_FILTER_DOUBLETS.out.summary.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_CLUSTERING:
    try {
      res_files = res_files.mix(ARCHR_CLUSTERING.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_EMBEDDING:
    try {
      res_files = res_files.mix(ARCHR_EMBEDDING.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MARKER_GENE:
    try {
      res_files = res_files.mix(ARCHR_MARKER_GENE.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_SCRNASEQ_UNCONSTRAINED:
    try {
      res_files = res_files.mix(ARCHR_SCRNASEQ_UNCONSTRAINED.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_SCRNASEQ_CONSTRAINED:
    try {
      res_files = res_files.mix(ARCHR_SCRNASEQ_CONSTRAINED.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_MARKER_PEAKS_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_GET_MARKER_PEAKS_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_MARKER_PEAKS_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_GET_MARKER_PEAKS_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_PAIRWISE_TEST_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_PAIRWISE_TEST_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_PAIRWISE_TEST_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_PAIRWISE_TEST_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MOTIF_ENRICHMENT_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_MOTIF_ENRICHMENT_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MOTIF_ENRICHMENT_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_MOTIF_ENRICHMENT_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MOTIF_DEVIATIONS_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_MOTIF_DEVIATIONS_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_MOTIF_DEVIATIONS_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_MOTIF_DEVIATIONS_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_COACCESSIBILITY_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_COACCESSIBILITY_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_COACCESSIBILITY_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_COACCESSIBILITY_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_TRAJECTORY_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_TRAJECTORY_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS:
    try {
      res_files = res_files.mix(ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2:
    try {
      res_files = res_files.mix(ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2.out.report.collect().ifEmpty([]))
    } catch (Exception ex) {}

  emit:
    res_files.collect()
    ARCHR_GET_CLUSTERING_TSV.out.res // Here if using collect(), only the first element will be used for split_bed module. For split bed
    ARCHR_GET_CLUSTERING_TSV.out.tsv // for split bam
}

// workflow.onComplete {
//     Completion.summary(workflow, params, log)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
