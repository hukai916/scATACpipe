// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CORRECT_BARCODE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode', publish_id:'') }
    container "hukai916/r_env_v4.4.1_amd64:0.1"

    input:
    tuple val(sample_name), path(barcode_fastq), path(valid_barcodes)

    output:
    tuple val(sample_name), val(chunk_name), path("tagfile_*.txt"), emit: sample_name_chunk_name_tagfile
    path "summary_*.txt", emit: corrected_barcode_summary

    script:
    chunk_name = barcode_fastq.name.split("\\.")[0..-3].join(".").split("barcode_").join() // get rid of suffix ".fastq.gz", then remove leading "barcode_" if there.

    """
    correct_barcode.R $options.args \
    --barcode_file=$barcode_fastq \
    --whitelist_file=$valid_barcodes \
    --path_output_fq=./

    """
}
