// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process SPLIT_BED {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'split_bed', publish_id:'') }

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
    tuple val(sample_name), path(fragment), path(tsv)

    output:
    path "split_*", emit: split_bed

    script:

    """
    tsv=($tsv)

    for (( i=0; i<\${#tsv[@]}; i++ )); do
      split_bed.py \${tsv[\$i]} $fragment
    done

    """
}
