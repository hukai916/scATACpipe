// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAGE_SAMPLE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'stage_sample', publish_id:'') }
    container "ubuntu:xenial"

    input:
    tuple val(sample_name), path(path_fastq_1), path(path_fastq_2), path(path_barcode)

    output:
    path path_fastq_1, emit: read1_fastq
    path path_barcode, emit: barcode_fastq
    path path_fastq_2, emit: read2_fastq
    val sample_name, emit: sample_name

    script:

    """
    touch test.txt

    """
}
