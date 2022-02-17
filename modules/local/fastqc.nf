// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTQC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'fastqc', publish_id:'') }
    container "hukai916/fastqc_0.11.9:0.1"

    input:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip
    val sample_name, emit: sample_name

    script:

    """
    fastqc $options.args --threads $task.cpus $read1_fastq $read2_fastq $barcode_fastq

    """
}
