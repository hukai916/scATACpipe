// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADD_BARCODE_TO_READ_CHUNKS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'add_barcode_to_read_chunks', publish_id:'') }
    container "hukai916/sinto_xenial:0.1"

    input:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)

    output:
    tuple val(sample_name), path("${read1_barcoded_fastq}"), path("${read2_barcoded_fastq}"), emit: reads_0
    // Below are for GET_WHITELIST_BARCODE
    val sample_name, emit: sample_name
    path barcode_fastq, emit: barcode_fastq

    script:
    read1_barcoded_fastq = read1_fastq.name.split("\\.")[0..-3].join(".") + ".barcoded.fastq.gz"
    read2_barcoded_fastq = read2_fastq.name.split("\\.")[0..-3].join(".") + ".barcoded.fastq.gz"

    """
    # use the first read length from fastq file to determine the length since -b is required by sinto.
    barcode_length=\$((zcat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    sinto barcode $options.args --barcode_fastq $barcode_fastq --read1 $read1_fastq --read2 $read2_fastq -b \$barcode_length

    """
}
