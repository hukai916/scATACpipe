/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
// NfcoreSchema pops compilation error.

// Validate input parameters
// WorkflowScatacpipe.initialise(params, log)

// TODO: add required input list here:
// Check input path parameters to see if they exist
def checkPathParamList = [ params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

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
include { GET_10XGENOMICS_FASTQ } from '../modules/local/get_10xgenomics_fastq'   addParams( options: modules['get_10xgenomics_fastq'] )
include { CELLRANGER_ATAC_COUNT } from '../modules/local/cellranger_atac_count'   addParams( options: modules['cellranger_atac_count'] )
include { GET_WHITELIST_BARCODE } from '../modules/local/get_whitelist_barcode'
include { CORRECT_BARCODE       } from '../modules/local/correct_barcode'         addParams( options: modules['correct_barcode'] )
include { CORRECT_BARCODE_PHENIQS } from '../modules/local/correct_barcode_pheniqs' addParams( options: modules['correct_barcode_pheniqs'] )
include { MATCH_READS           } from '../modules/local/match_reads'             addParams( options: modules['match_reads'] )
include { MATCH_READS_TRIMMED   } from '../modules/local/match_reads_trimmed'     addParams( options: modules['match_reads_trimmed'] )
include { FASTQC                } from '../modules/local/fastqc'                  addParams( options: modules['fastqc'] )
include { ADD_BARCODE_TO_READS       } from '../modules/local/add_barcode_to_reads'    addParams( options: modules['add_barcode_to_reads'] )
include { CUTADAPT         } from '../modules/local/cutadapt'    addParams( options: modules['cutadapt'] )
include { DOWNLOAD_FROM_UCSC } from '../modules/local/download_from_ucsc'    addParams( options: modules['download_from_ucsc'] )
include { DOWNLOAD_FROM_ENSEMBL } from '../modules/local/download_from_ensembl'    addParams( options: modules['download_from_ensembl'] )
include { BWA_INDEX        } from '../modules/local/bwa_index'    addParams( options: modules['bwa_index'] )
include { BWA_MAP          } from '../modules/local/bwa_map'    addParams( options: modules['bwa_map'] )
include { PREP_GENOME } from '../modules/local/prep_genome'
include { BUILD_GENE_ANNOTATION } from '../modules/local/build_gene_annotation' addParams( options: modules['build_gene_annotation'] )
include { BUILD_GENOME_ANNOTATION } from '../modules/local/build_genome_annotation' addParams( options: modules['build_genome_annotation'] )
include { MINIMAP2_INDEX   } from '../modules/local/minimap2_index'    addParams( options: modules['minimap2_index'] )
include { MINIMAP2_MAP     } from '../modules/local/minimap2_map'    addParams( options: modules['minimap2_map'] )
include { BAM_FILTER       } from '../modules/local/bam_filter'    addParams( options: modules['bam_filter'] )
include { REMOVE_DUPLICATE } from '../modules/local/remove_duplicate'    addParams( options: modules['remove_duplicate'] )
include { QUALIMAP         } from '../modules/local/qualimap'    addParams( options: modules['qualimap'] )
include { GET_FRAGMENTS    } from '../modules/local/get_fragments'    addParams( options: modules['get_fragments'] )
include { DOWNLOAD_FROM_UCSC_GTF } from '../modules/local/download_from_ucsc_gtf'    addParams( options: modules['download_from_ucsc_gtf'] )
include { DOWNLOAD_FROM_ENSEMBL_GTF } from '../modules/local/download_from_ensembl_gtf'    addParams( options: modules['download_from_ensembl_gtf'] )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow PREPROCESS_DEFAULT {
  take:
    reads

  main:
    // Examine if all required parameters supplied:
    if (!(params.mapper == "bwa") && !(params.mapper == "minimap2")) {
      log.error "Must supply --mapper [bwa | minimap2]!"
      exit 1, "EXIT!"
    } else if (!(params.ref_bwa_index) && !(params.ref_minimap2_index) && !(params.ref_fasta) && !(params.ref_fasta_ensembl) && !(params.ref_fasta_ucsc)) {
      msg = "Must supply one from below: \n" + "  --ref_bwa_index\n" + "  --ref_minimap2_index\n" + "  --ref_fasta\n" + "  --ref_fasta_ensembl\n" + "  --ref_fasta_ucsc\n"
      log.error msg
      exit 1, "EXIT!"
    }
    // Above is redundant to WorkflowMain::initialise()

    // log.info "INFO(2): --preprocess: default"
    GET_10XGENOMICS_FASTQ (reads)
    // module: fastQC
    FASTQC (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq,
    GET_10XGENOMICS_FASTQ.out.barcode_fastq)

    // Module: barcode correction (optional) and add barcode: correct barcode fastq given whitelist and barcode fastq file
    if (!(params.barcode_whitelist)) {
      log.info "NOTICE: --barcode_whitelist: not supplied, skip barcode correction!"
      ADD_BARCODE_TO_READS (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
    } else {
      // Module: determine the right whitelist barcode
      GET_WHITELIST_BARCODE (GET_10XGENOMICS_FASTQ.out.sample_name, GET_10XGENOMICS_FASTQ.out.barcode_fastq, Channel.fromPath(params.barcode_whitelist).collect(), GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)

      // Allow users to choose from barcode_correction.R or Pheniqs:
      if (params.barcode_correction == "naive") {
        CORRECT_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, GET_WHITELIST_BARCODE.out.whitelist_barcode, GET_WHITELIST_BARCODE.out.read1_fastq, GET_WHITELIST_BARCODE.out.read2_fastq)
      // MATCH_READS (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.corrected_barcode, GET_10XGENOMICS_FASTQ.out.read1_fastq, GET_10XGENOMICS_FASTQ.out.read2_fastq)
      // Note that the above might be problematic, since MATCH_READS would take inputs from two channels, the instance of samples may not match.
        MATCH_READS (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.corrected_barcode, CORRECT_BARCODE.out.read1_fastq, CORRECT_BARCODE.out.read2_fastq)

        ADD_BARCODE_TO_READS (MATCH_READS.out.sample_name, MATCH_READS.out.barcode1_fastq, MATCH_READS.out.barcode2_fastq, MATCH_READS.out.read1_fastq, MATCH_READS.out.read2_fastq)
      } else {
        // use pheniqs:
        CORRECT_BARCODE_PHENIQS (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, GET_WHITELIST_BARCODE.out.whitelist_barcode, GET_WHITELIST_BARCODE.out.read1_fastq, GET_WHITELIST_BARCODE.out.read2_fastq)

        MATCH_READS (CORRECT_BARCODE_PHENIQS.out.sample_name, CORRECT_BARCODE_PHENIQS.out.corrected_barcode, CORRECT_BARCODE_PHENIQS.out.read1_fastq, CORRECT_BARCODE_PHENIQS.out.read2_fastq)
      }
    }

    // module: trimming off adapter
    if ((params.barcode_whitelist) && (params.barcode_correction == "pheniqs")) {
      CUTADAPT (MATCH_READS.out.sample_name, MATCH_READS.out.read1_fastq, MATCH_READS.out.read2_fastq, params.read1_adapter, params.read2_adapter)
    } else {
      CUTADAPT (ADD_BARCODE_TO_READS.out.sample_name, ADD_BARCODE_TO_READS.out.read1_fastq, ADD_BARCODE_TO_READS.out.read2_fastq, params.read1_adapter, params.read2_adapter)
    }

    // module: MATCH_READS_TRIMMED: in case user choose to trim based on quality and read pair gets unbalanced.
    MATCH_READS_TRIMMED (CUTADAPT.out.sample_name, CUTADAPT.out.trimed_read1_fastq, CUTADAPT.out.trimed_read2_fastq)

    // module: mapping with bwa or minimap2: mark duplicate
    // bwa or minimap2
    if (params.mapper == 'bwa') {
      log.info "INFO: --mapper: bwa"
      if (params.ref_bwa_index) {
        BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, params.ref_bwa_index)
      } else if (params.ref_fasta) {
        log.info "INFO: --ref_fasta provided, use it for building bwa index."
        // module : prep_genome
        PREP_GENOME (params.ref_fasta, "custom_genome")
        // module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
        // module : bwa_map
        BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, BWA_INDEX.out.bwa_index_folder.collect())
      } else if (params.ref_fasta_ensembl) {
        log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // module : download_from_ensembl
        DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // module: prep_genome
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
      } else if (params.ref_fasta_ucsc) {
        log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build bwa index, and map with bwa ..."
        // module : download_from_ucsc
        DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // module : prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_name)
        // module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
        // module : bwa_map
        BWA_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, BWA_INDEX.out.bwa_index_folder.collect())
      } else {
        exit 1, 'Parameter --ref_fasta_ensembl/--ref_fasta_ucsc: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
      }
    } else if (params.mapper == "minimap2") {
      log.info "INFO: --mapper: minimap2"
      if (params.ref_minimap2_index) {
        // use user provided bwa index for mapping
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, params.ref_minimap2_index)
      } else if (params.ref_fasta) {
        log.info "INFO: --ref_fasta provided, use it to build minimap2 index."
        // module : prep_genome
        PREP_GENOME (params.ref_fasta, "custom_genome")
        // module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, MINIMAP2_INDEX.out.minimap2_index.collect())
      } else if (params.ref_fasta_ensembl) {
        log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // module : download_from_ensembl
        DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // module: PREP_GENOME
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, MINIMAP2_INDEX.out.minimap2_index.collect())
      } else if (params.ref_fasta_ucsc) {
        log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // module : download_from_ucsc
        DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // module : prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_gtf)
        // module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // module : minimap2_map
        MINIMAP2_MAP (MATCH_READS_TRIMMED.out.sample_name, MATCH_READS_TRIMMED.out.read1_fastq, MATCH_READS_TRIMMED.out.read2_fastq, MINIMAP2_INDEX.out.minimap2_index.collect())
      } else {
        exit 1, 'Parameter --ref_fasta_ucsc/--ref_fasta_ensembl: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
      }
    } else {
      log.error "--mapper must be supplied!"
      exit 1, "EXIT!"
    }

    // module: filter out poorly mapped reads
    if (params.mapper == 'bwa') {
      BAM_FILTER (BWA_MAP.out.sample_name, BWA_MAP.out.bam, params.filter)
    } else if (params.mapper == "minimap2") {
        BAM_FILTER (MINIMAP2_MAP.out.sample_name, MINIMAP2_MAP.out.bam, params.filter)
    }

    // module: remove duplicates based on cell barcode, start, end
    REMOVE_DUPLICATE(BAM_FILTER.out.sample_name, BAM_FILTER.out.bam)

    // module: bamqc with qualimap for raw bam files
    // TODO: Run Qualimap on the final filtered deduplicated bam file.
    QUALIMAP (REMOVE_DUPLICATE.out.sample_name, REMOVE_DUPLICATE.out.bam)

    // module: generate fragment file with sinto
    // use raw bam file since ArchR may take advantage of the duplication info.
    GET_FRAGMENTS (BAM_FILTER.out.sample_name, BAM_FILTER.out.bam)

    // Collect all output results for MultiQC report:
    res_files = Channel.empty()
    // res_files = res_files.mix(Channel.from(ch_multiqc_config))
    // res_files = res_files.mix(Channel.from(ch_multiqc_custom_config).collect().ifEmpty([]))

    // Use try-catch since if certain module is not run, module.out becomes undefined.
    // FASTQC module:
    try {
      res_files = res_files.mix(FASTQC.out.zip.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // CORRECT_BARCODE module:
    try {
      res_files = res_files.mix(CORRECT_BARCODE.out.corrected_barcode_summary.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // CORRECT_BARCODE_PHENIQS module:
    try {
      res_files = res_files.mix(CORRECT_BARCODE_PHENIQS.out.corrected_barcode_summary.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // REMOVE_DUPLICATE module:
    try {
      res_files = res_files.mix(REMOVE_DUPLICATE.out.remove_duplicate_summary.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // CUTADAPT module:
    try {
      res_files = res_files.mix(CUTADAPT.out.log.collect().ifEmpty([]))
    } catch (Exception ex) {}
    // QUALIMAP module:
    try {
      res_files = res_files.mix(QUALIMAP.out.bamqc.collect().ifEmpty([]))
    } catch (Exception ex) {}

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
      // prep_gtf          = Channel.value("run")
    } catch (Exception ex) {}
    try {
      prep_gtf_file = PREP_GTF.out.gtf
    } catch (Exception ex) { }

  emit:
    res_files // out[0]: res folders for MultiQC report
    GET_FRAGMENTS.out.fragments // out[1]: for split bed
    GET_FRAGMENTS.out.ch_fragment // out[2]: fragment ch for ArchR
    REMOVE_DUPLICATE.out.sample_name // out[3]: for split bam
    REMOVE_DUPLICATE.out.bam // out[4]: for split bam
    prep_genome_name         // out[5]: for DOWNSTREAM_ARCHR
    prep_genome_fasta        // out[6]: for DOWNSTREAM_ARCHR
    prep_gtf_genome          // out[7]: for DOWNSTREAM_ARCHR
    prep_gtf_file            // out[8]: for DOWNSTREAM_ARCHR
}

// workflow.onComplete {
//     Completion.summary(workflow, params, log)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
