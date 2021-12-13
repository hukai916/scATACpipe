// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CORRECT_BARCODE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode', publish_id:'') }
    container "hukai916/r_sc_atac:0.1"

    input:
    // val sample_name
    // path read1_fastq
    // path read2_fastq
    // path barcode_fastq
    // path barcode_whitelist
    val sample_name
    path barcode_fastq
    path valid_barcodes


    output:
    val sample_name, emit: sample_name
    path "tagfile_*.txt", emit: tagfile
    path "summary_*.txt", emit: corrected_barcode_summary

    // tuple val(sample_name), path(read1_fastq), path(read2_fastq), path("barcode_*"), emit: reads
    // path "summary_*.txt", emit: corrected_barcode_summary

    script:

    """
    correct_barcode.R $options.args \
    --barcode_file=$barcode_fastq \
    --whitelist_file=$valid_barcodes \
    --path_output_fq=./

    """
}
