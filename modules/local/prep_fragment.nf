// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process PREP_FRAGMENT {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'prep_fragment', publish_id:'') }

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
    tuple val(sample_name), path(fragment)
    path gtf

    output:
    tuple val(sample_name), path("final.fragment.bed.gz"), emit: fragment

    """
    if [[ $fragment == *.gz ]]; then
      gunzip -c $fragment > fragment.bed
    else
      if [[ $fragment != fragment.bed ]]; then
        mv $fragment fragment.bed
      fi
    fi

    if [[ $gtf == *.gz ]]; then
      gunzip -c $gtf > annotation.gtf
    else
      if [[ $gtf != annotation.gtf ]]; then
        mv $gtf annotation.gtf
      fi
    fi

    # add 'chr' prefix to fragment bed to match PREP_GTF
    cat fragment.bed | awk 'BEGIN{FS=OFS="\\t"}/^#/{print; next} NF>=3{\$1 = gensub(/(chr)?(.+).*/, "chr\\\\2", "g", \$1); print}' > chrPrefixed.fragment.bed

    # make sure bed col1 is a subset of gtf col1:
    extract_bed.py chrPrefixed.fragment.bed annotation.gtf | bgzip > final.fragment.bed.gz

    """
}
