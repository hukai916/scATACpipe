// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process MINIMAP2_MAP {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'minimap2_map', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/minimap2_xenial:0.1"

    // cache false

    input:
    val sample_name
    path read1_fastq
    path read2_fastq
    path minimap2_index_file

    output:
    val sample_name, emit: sample_name
    path "*.sorted.bam", emit: bam

    script:

    """
    minimap2 -t $task.cpus $options.args -a $minimap2_index_file $read1_fastq $read2_fastq | samtools sort -@ $task.cpus -O bam -o ${sample_name}.sorted.bam
    # note that -ax sr pops error for test dataset.
    # minimap2 is not good for short reads alignments.
    """
}
