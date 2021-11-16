// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process GET_PRIMARY_GENOME {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_primary_genome', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/bwa_xenial:0.1"
    // cache false

    input:
    path genome_fasta

    output:
    path "primary_genome.fa.gz", emit: genome_fasta

    """
    gunzip -c $genome_fasta > genome.fa
    samtools faidx genome.fa

    cat genome.fa.fai | cut -f 1 | grep -v "_alt\$|_fix\$" | xargs -n 1 -I {} samtools faidx genome.fa {} >> primary_genome.fa
    gzip primary_genome.fa

    rm genome.fa genome.fa.fai
    """
}
