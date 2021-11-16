// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process REMOVE_DUPLICATE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'remove_duplicate', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/pysam_xenial:0.1"

    // cache false

    input:
    val sample_name
    path bam

    output:
    val sample_name, emit: sample_name
    path "rm_dup_*.sorted.bam", emit: bam
    path "summary_rm_dup_*.txt", emit: remove_duplicate_summary

    script:

    """
    # remove PCR duplicates based on cell barcode, start, end:
    remove_duplicate.py --inbam $bam --outbam rm_dup_${sample_name}.bam $options.args

    # sort output bam:
    samtools sort rm_dup_${sample_name}.bam -o rm_dup_${sample_name}.sorted.bam

    """

}
