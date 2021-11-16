// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process FASTQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'fastqc', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/fastqc_0.11.9:0.1"
    // cache false

    input:
    val sample_name
    path read1_fastq
    path read2_fastq
    path barcode_fastq

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip
    val sample_name, emit: sample_name

    script:

    """
    fastqc $options.args --threads $task.cpus $read1_fastq $read2_fastq $barcode_fastq

    """
}
