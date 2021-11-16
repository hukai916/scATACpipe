// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process FIX_UCSC_GTF {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'fix_ucsc_gtf', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/miniconda3_bio:0.1"

    // cache false

    input:
    path gtf
    path genome

    output:
    path "fixed.*.gz", emit: gtf

    script:

    """
    gunzip -c $gtf > annotation.gtf
    # sort the gtf in case some records are not in order
    sort_gtf.py annotation.gtf > sorted.annotation.gtf

    # make sure gtf config is a subset of genome.fa config, required for cellranger_index step
    extract_gtf.py $genome sorted.annotation.gtf > fixed.annotation.gtf

    gzip fixed.annotation.gtf

    """
}
