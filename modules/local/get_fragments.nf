// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process GET_FRAGMENTS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'fragments', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/sinto_xenial:0.2"

    // cache false

    input:
    val sample_name
    path bam

    output:
    val sample_name, emit: sample_name
    path "fragments.sort.bed.gz", emit: fragments
    tuple val(sample_name), path("fragments.sort.bed.gz"), emit: ch_fragment

    script:

    """
    # first index the bam file
    samtools index $options.args $bam

    # then, generate the fragments file
    sinto fragments $options.args --nproc $task.cpus --bam $bam -f fragments.bed --barcode_regex "[^:]*"
    # sort and bzip the fragment file
    sort -k 1,1 -k2,2n fragments.bed > fragments.sort.bed
    bgzip fragments.sort.bed

    """
}
