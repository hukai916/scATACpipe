// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
include { get_genome_ucsc } from './genome_ucsc'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process DOWNLOAD_TEST {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_test', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/miniconda_bio:0.1"

    // cache false

    input:
    val test

    output:
    // path "text.txt"

    script:
    dict_genome_ucsc = get_genome_ucsc()
    dict_genome_ucsc.each { println "Supported UCSC genome: $it.key" }

    """
    """
}
