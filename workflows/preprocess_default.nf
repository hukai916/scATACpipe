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
include { SPLIT_FASTQ           } from '../modules/local/split_fastq'
include { MATCH_CHUNK           } from '../modules/local/match_chunk'
include { GET_WHITELIST_BARCODE } from '../modules/local/get_whitelist_barcode'
include { GET_VALID_BARCODE } from '../modules/local/get_valid_barcode'
include { CORRECT_BARCODE       } from '../modules/local/correct_barcode'         addParams( options: modules['correct_barcode'] )
include { CORRECT_BARCODE_PHENIQS } from '../modules/local/correct_barcode_pheniqs' addParams( options: modules['correct_barcode_pheniqs'] )
include { MATCH_READS           } from '../modules/local/match_reads'             addParams( options: modules['match_reads'] )
include { MATCH_READS_TRIMMED   } from '../modules/local/match_reads_trimmed'     addParams( options: modules['match_reads_trimmed'] )
include { FASTQC                } from '../modules/local/fastqc'                  addParams( options: modules['fastqc'] )
include { ADD_BARCODE_TO_READS       } from '../modules/local/add_barcode_to_reads'
include { ADD_BARCODE_TO_READS_2       } from '../modules/local/add_barcode_to_reads_2'
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
include { FILTER_BAM       } from '../modules/local/filter_bam'    addParams( options: modules['filter_bam'] )
// include { PREP_BAM         } from '../modules/local/prep_bam'    addParams( options: modules['prep_bam'] )
include { TAG_BAM          } from '../modules/local/tag_bam'
include { DEDUP_BAM; DEDUP_BAM as DEDUP_BAM2 } from '../modules/local/dedup_bam'
// include { ADD_BARCODE_TO_TAG } from '../modules/local/add_barcode_to_tag'
include { QUALIMAP         } from '../modules/local/qualimap'    addParams( options: modules['qualimap'] )
include { GET_FRAGMENTS    } from '../modules/local/get_fragments'    addParams( options: modules['get_fragments'] )
include { DOWNLOAD_FROM_UCSC_GTF } from '../modules/local/download_from_ucsc_gtf'    addParams( options: modules['download_from_ucsc_gtf'] )
include { DOWNLOAD_FROM_ENSEMBL_GTF } from '../modules/local/download_from_ensembl_gtf'    addParams( options: modules['download_from_ensembl_gtf'] )
include { COMBINE_FRAGMENTS } from '../modules/local/combine_fragments'
include { COMBINE_BAM; COMBINE_BAM as COMBINE_BAM2 } from '../modules/local/combine_bam'


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow PREPROCESS_DEFAULT {
  take:
    reads
    sample_count

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

    // .out.reads: tuple val(sample_name), path(read1), path(read2), path(barcode)
    // .out.reads_0: tuple val(sample_name), path(read1), path(read2)
    // .out.reads_2: tuple val(sample_name), path(read1), path(read2), path(barcode1), path(barcode2)

    // log.info "INFO(2): --preprocess: default"

    // Module: fastQC
    FASTQC (reads)

    if (!params.split_fastq) {
      ADD_BARCODE_TO_READS (reads)
    } else {
      // Module: split read into 20M chunks
      SPLIT_FASTQ (reads, sample_count)
      // Below is to ensure that sample_name, R1, R2, and barcode matches
      read1_chunk    = SPLIT_FASTQ.out.read1_fastq.collect().flatten()
      read2_chunks   = SPLIT_FASTQ.out.read2_fastq.collect()
      barcode_chunks = SPLIT_FASTQ.out.barcode_fastq.collect()
      MATCH_CHUNK (read1_chunk, read2_chunks, barcode_chunks)
      sample_name    = MATCH_CHUNK.out.sample_name.unique() // note .collect().unique() won't do the job probably because collect() return a single list, and unique() must work on channel.

      // Module: add barcode to reads
      ADD_BARCODE_TO_READS (MATCH_CHUNK.out.chunk)
    }

    // Module: trim off adapter
    CUTADAPT (ADD_BARCODE_TO_READS.out.reads_0, params.read1_adapter, params.read2_adapter)

    // Module: mapping with bwa or minimap2: mark duplicate
    // bwa or minimap2
    if (params.mapper == 'bwa') {
      log.info "INFO: --mapper: bwa"
      if (params.ref_bwa_index) {
        BWA_MAP (CUTADAPT.out.reads_0, params.ref_bwa_index)
      } else if (params.ref_fasta) {
        log.info "INFO: --ref_fasta provided, use it for building bwa index."
        // Module : prep_genome
        PREP_GENOME (params.ref_fasta, "custom_genome")
        // Module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
        // Module : bwa_map
        BWA_MAP (CUTADAPT.out.reads_0, BWA_INDEX.out.bwa_index_folder.collect())
      } else if (params.ref_fasta_ensembl) {
        log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // Module : download_from_ensembl
        DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // Module: prep_genome
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // Module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
      } else if (params.ref_fasta_ucsc) {
        log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build bwa index, and map with bwa ..."
        // Module : download_from_ucsc
        DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // Module : prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_name)
        // Module : bwa_index
        BWA_INDEX (PREP_GENOME.out.genome_fasta)
        // Module : bwa_map
        BWA_MAP (CUTADAPT.out.reads_0, BWA_INDEX.out.bwa_index_folder.collect())
      } else {
        exit 1, 'Parameter --ref_fasta_ensembl/--ref_fasta_ucsc: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
      }
    } else if (params.mapper == "minimap2") {
      log.info "INFO: --mapper: minimap2"
      if (params.ref_minimap2_index) {
        // use user provided bwa index for mapping
        // Module : minimap2_map
        MINIMAP2_MAP (CUTADAPT.out.reads_0, params.ref_minimap2_index)
      } else if (params.ref_fasta) {
        log.info "INFO: --ref_fasta provided, use it to build minimap2 index."
        // Module : prep_genome
        PREP_GENOME (params.ref_fasta, "custom_genome")
        // Module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // Module : minimap2_map
        MINIMAP2_MAP (CUTADAPT.out.reads_0, MINIMAP2_INDEX.out.minimap2_index.collect())
      } else if (params.ref_fasta_ensembl) {
        log.info "INFO: --ref_fasta_ensembl provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // Module : download_from_ensembl
        DOWNLOAD_FROM_ENSEMBL (params.ref_fasta_ensembl, Channel.fromPath('assets/genome_ensembl.json'))
        // Module: PREP_GENOME
        PREP_GENOME (DOWNLOAD_FROM_ENSEMBL.out.genome_fasta, DOWNLOAD_FROM_ENSEMBL.out.genome_name)
        // Module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // Module : minimap2_map
        MINIMAP2_MAP (CUTADAPT.out.reads_0, MINIMAP2_INDEX.out.minimap2_index.collect())
      } else if (params.ref_fasta_ucsc) {
        log.info "INFO: --ref_fasta_ucsc provided, will download genome, and then build minimap2 index, and map with minimap2 ..."
        // Module : download_from_ucsc
        DOWNLOAD_FROM_UCSC (params.ref_fasta_ucsc, Channel.fromPath('assets/genome_ucsc.json'))
        // Module : prep_genome
        PREP_GENOME (DOWNLOAD_FROM_UCSC.out.genome_fasta, DOWNLOAD_FROM_UCSC.out.genome_gtf)
        // Module : bwa_index
        MINIMAP2_INDEX (PREP_GENOME.out.genome_fasta)
        // Module : minimap2_map
        MINIMAP2_MAP (CUTADAPT.out.reads_0, MINIMAP2_INDEX.out.minimap2_index.collect())
      } else {
        exit 1, 'Parameter --ref_fasta_ucsc/--ref_fasta_ensembl: pls supply a genome name, like hg19, mm10 (if ucsc), or homo_sapiens, mus_musculus (if ensembl)!'
      }
    } else {
      log.error "--mapper must be supplied!"
      exit 1, "EXIT!"
    }

    // Module: filter out poorly mapped reads
    if (params.mapper == 'bwa') {
      FILTER_BAM (BWA_MAP.out.sample_name, BWA_MAP.out.bam, params.filter)
    } else if (params.mapper == "minimap2") {
      FILTER_BAM (MINIMAP2_MAP.out.sample_name, MINIMAP2_MAP.out.bam, params.filter)
    }

    if (!params.split_fastq) {
      // Module: dedup bam by barcode seq added in the front
      DEDUP_BAM (FILTER_BAM.out.sample_name, FILTER_BAM.out.bam, "N/A")
    } else {
      // Module: combine bam
      COMBINE_BAM (FILTER_BAM.out.sample_name.unique(), FILTER_BAM.out.bam.collect())
      // Module: dedup bam by barcode seq added in the front
      DEDUP_BAM (COMBINE_BAM.out.sample_name, COMBINE_BAM.out.bam, "N/A")
    }

    if (!params.barcode_correction) {
      // Module: get fragment file
      GET_FRAGMENTS (DEDUP_BAM.out.sample_name, DEDUP_BAM.out.bam)
    } else {
      // Note, we first obtain valid barcode on a per sample basis no matter split_fastq or not:
      // Module: get_whitelist_barcode
      if (!params.whitelist_barcode) {
        use_whitelist = "false"
        GET_WHITELIST_BARCODE (reads, Channel.fromPath('assets/whitelist_barcodes').first())
      } else {
        use_whitelist = "true"
        reads.view()
        log.info "test here"
        GET_WHITELIST_BARCODE (reads, Channel.fromPath(params.whitelist_barcode).first())
      }
      // Module: get_valid_barcode
      GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.barcode_whitelist.join(DEDUP_BAM.out.bam), use_whitelist) // sample_name, whitelist_barcode

      // Module: get valid barcode
      // if (!params.split_fastq) {
      //   if (!params.whitelist_barcode) {
      //     use_whitelist = "false"
      //     // per sample, need only one whitelist barcode:
      //     GET_WHITELIST_BARCODE (reads, Channel.fromPath('assets/whitelist_barcodes').first())
      //     // per sample, need only one valid barcode:
      //     //sample_name, whitelist_barcode, raw_barcode_fastq, dedup_bam
      //     GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.whitelist_barcode, use_whitelist) // sample_name, whitelist_barcode
      //
      //     // A cross operator is needed here: DEDUP_BAM::sample_name with GET_WHITELIST_BARCODE::sample_name
      //     // Need: sample_name, raw barcodes, dedup_bam, whitelist; Output: barcode_fastq for that read/chunk, sample_name, valid_counts
      //     GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, )
      //
      //     GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, DEDUP_BAM.out.bam.collect(), Channel.fromPath("assets/file_token.txt").first())
      //
      //     GET_WHITELIST_BARCODE (sample_name, barcode_chunks, Channel.fromPath('assets/whitelist_barcodes').first())
      //     GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, DEDUP_BAM.out.bam.collect(), Channel.fromPath("assets/file_token.txt").first())
      //   } else {
      //     GET_WHITELIST_BARCODE (sample_name, barcode_chunks, Channel.fromPath(params.whitelist_barcode).first())
      //     GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, DEDUP_BAM.out.bam.collect(), GET_WHITELIST_BARCODE.out.whitelist_barcode)
      //   }
      // } else {
      //   if (!params.whitelist_barcode) {
      //     GET_WHITELIST_BARCODE (sample_name, barcode_chunks, Channel.fromPath('assets/whitelist_barcodes').first())
      //     GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, DEDUP_BAM.out.bam.collect(), Channel.fromPath("assets/file_token.txt").first())
      //   } else {
      //     GET_WHITELIST_BARCODE (sample_name, barcode_chunks, Channel.fromPath(params.whitelist_barcode).first())
      //     GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, DEDUP_BAM.out.bam.collect(), GET_WHITELIST_BARCODE.out.whitelist_barcode)
      //   }
      // }
      // Module: get valid barcode
      // if (!params.whitelist_barcode) {
      //   GET_WHITELIST_BARCODE (sample_name, barcode_chunks, Channel.fromPath('assets/whitelist_barcodes').first())
      //   GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, DEDUP_BAM.out.bam.collect(), Channel.fromPath("assets/file_token.txt").first())
      // } else {
      //   GET_WHITELIST_BARCODE (sample_name, barcode_chunks, Channel.fromPath(params.whitelist_barcode).first())
      //   GET_VALID_BARCODE (GET_WHITELIST_BARCODE.out.sample_name, GET_WHITELIST_BARCODE.out.barcode_fastq, DEDUP_BAM.out.bam.collect(), GET_WHITELIST_BARCODE.out.whitelist_barcode)
      // }

      // Prepare for CORRECT_BARCODE input:
      if (!params.split_fastq) {
        left = reads
      } else {
        left = MATCH_CHUNK.out.chunk
      }
      left
        .combine(GET_VALID_BARCODE.out.valid_barcodes_and_counts)
        .filter({ it[0] == it[4] })
        .set(ch_correct_barcode_input)
      ch_correct_barcode_input
        .map({ it -> [ it[0], it[3], it[5]] })
        .set(ch_correct_barcode_pheniqs_input) // sample_name, barcode_fastq, valid_barcode_couts_fastq
      ch_correct_barcode_input
        .map({ it -> [ it[0], it[3], it[4]] })
        .set(ch_correct_barcode_naive_input) // sample_name, barcode_fastq, valid_barcodes

      // Modules: correct barcode
      if (params.barcode_correction == "pheniqs") {
      // Module: pheniqs
        CORRECT_BARCODE_PHENIQS (ch_correct_barcode_pheniqs_input)
      // Module: tag_bam
        CORRECT_BARCODE_PHENIQS
          .out
          .sample_name_chunk_name_tagfile
          .join(FILTER_BAM.out.sample_name_chunk_name_bam, by: [0, 1])
          .set(ch_tag_bam_input) // sample_name, file_name, tagfile, bam
      // Module: add CB tag to BAM containg corrected barcodes
        TAG_BAM (ch_tag_bam_input)
      } else if (params.barcode_correction == "naive") {
        // Module: naive barcode correction
        // CORRECT_BARCODE (sample_name, read1_chunk, read2_chunk, barcode_chunk, GET_VALID_BARCODE.out.valid_barcode.collect())

        // NEEDs update:
        CORRECT_BARCODE (GET_VALID_BARCODE.out.sample_name, GET_VALID_BARCODE.out.barcode_fastq, GET_VALID_BARCODE.out.valid_barcodes)

        // Module: add CB tag to BAM containing corrected barcodes
        // TAG_BAM (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.tagfile, DEDUP_BAM.out.bam.collect())
        TAG_BAM (CORRECT_BARCODE.out.sample_name, CORRECT_BARCODE.out.tagfile, COMBINE_BAM.out.bam.collect())
      }
      // Module: dedup bam again using "CB" tag
      if (!params.split_fastq) {
        DEDUP_BAM2 (TAG_BAM.out.sample_name, TAG_BAM.out.bam, "CB")
      } else {
        COMBINE_BAM2 (TAG_BAM.out.sample_name.unique(), TAG_BAM.out.bam.collect())
        DEDUP_BAM2 (COMBINE_BAM2.out.sample_name, COMBINE_BAM2.out.bam, "CB")
      }
      // Module: get fragments
      GET_FRAGMENTS (DEDUP_BAM2.out.sample_name, DEDUP_BAM2.out.bam)
    }

    // Module: barcode correction (optional) and add barcode: correct barcode fastq given whitelist and barcode fastq file
    // Module: remove duplicates based on cell barcode, start, end
    // REMOVE_DUPLICATE(COMBINE_BAM.out.sample_name, COMBINE_BAM.out.bam)
    // DISCUSS: bamqc with qualimap for raw bam files
    // QUALIMAP (REMOVE_DUPLICATE.out.sample_name, REMOVE_DUPLICATE.out.bam)

    // Module: generate fragment file with sinto
    // use raw bam file since ArchR may take advantage of the duplication info.
    // GET_FRAGMENTS (REMOVE_DUPLICATE.out.sample_name, REMOVE_DUPLICATE.out.bam)
    // GET_FRAGMENTS (DEDUP_BAM2.out.sample_name, DEDUP_BAM2.out.bam)

    // Module: run Qualimap on the final filtered, deduplicated, combined, and sorted bam file.
    QUALIMAP (DEDUP_BAM2.out.sample_name, DEDUP_BAM2.out.bam)

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
    // COMBINE_FRAGMENTS.out.fragments // out[1]: for split bed
    // COMBINE_FRAGMENTS.out.ch_fragment // out[2]: fragment ch for ArchR
    GET_FRAGMENTS.out.fragments // out[1]: for split bed
    GET_FRAGMENTS.out.ch_fragment // out[2]: fragment ch for ArchR
    COMBINE_BAM.out.sample_name // out[3]: for split bam
    COMBINE_BAM.out.bam // out[4]: for split bam
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
