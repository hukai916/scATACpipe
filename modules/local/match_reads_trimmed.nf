// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process MATCH_READS_TRIMMED {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'match_reads_trimmed', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/seqkit_0.16.1:0.1"

    // cache false

    input:
    val sample_name
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name
    path "R1/*.fastq.gz", emit: read1_fastq
    path "R2/*.fastq.gz", emit: read2_fastq

    script:

    """
    seqkit pair $options.args -1 $read1_fastq -2 $read2_fastq -O R1
    mkdir R2
    mv R1/$read2_fastq R2/

    """
}
