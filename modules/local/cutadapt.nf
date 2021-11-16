// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process CUTADAPT {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'cutadapt', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/cutadapt_xenial:0.1"

    // cache false

    input:
    val sample_name
    path read1_fastq
    path read2_fastq
    val read1_adapter
    val read2_adapter

    output:
    val sample_name, emit: sample_name
    path "R1/trimmed*", emit: trimed_read1_fastq
    path "R2/trimmed*", emit: trimed_read2_fastq
    path "log_cutadapt_*.txt", emit: log

    script:
    read1_name = read1_fastq.getName()
    read2_name = read2_fastq.getName()

    """
    mkdir R1
    cutadapt $options.args -a $read1_adapter -o R1/trimmed_$read1_name $read1_fastq
    mkdir R2
    cutadapt $options.args -a $read2_adapter -o R2/trimmed_$read2_name $read2_fastq
    cp .command.log log_cutadapt_${sample_name}.txt

    """
}
