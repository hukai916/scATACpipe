// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process CELLRANGER_ATAC_COUNT {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'cellranger_count', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/cellranger_atat_2.0.0:0.1"

    // cache false

    input:
    val sample_name
    path fastq_folder
    path reference
    // val jobmode

    output:
    val sample_name, emit: sample_name
    path "cellranger_atac_count_*/outs/fragments.tsv.gz", emit: fragments
    tuple val(sample_name), path("cellranger_atac_count_*/outs/fragments.tsv.gz"), emit: ch_fragment
    path "cellranger_atac_count_*", emit: cellranger_atac_count
    path "cellranger_atac_count_*/outs/*_possorted_bam.bam", emit: bam

    script:
    // def avail_mem = task.memory ? "${ (task.memory.toBytes().intdiv(1073741824).intdiv(task.cpus) * 0.9).toInteger() }" : ''
    def avail_mem = task.memory ? "${ (task.memory.toBytes().intdiv(1073741824) * 0.9).toInteger() }" : ''

    """
    # the fastq file name must not contain special characters other than dash, underscore, digit; dot is not allowed
    # the --id must not contain dot either:

    fastq_folder=\$( echo $fastq_folder | tr '.' '_' ) # just in case

    cellranger-atac count $options.args \
    --id cellranger_atac_count_\$fastq_folder \
    --fastqs $fastq_folder \
    --reference $reference \
    --localcores $task.cpus \
    --localmem $avail_mem

    # rename the output bam file for split_bam module:
    mv cellranger_atac_count_\${fastq_folder}/outs/possorted_bam.bam cellranger_atac_count_\${fastq_folder}/outs/${sample_name}_possorted_bam.bam
    mv cellranger_atac_count_\${fastq_folder}/outs/possorted_bam.bam.bai cellranger_atac_count_\${fastq_folder}/outs/${sample_name}_possorted_bam.bam.bai

    """
}
