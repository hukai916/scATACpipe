// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_FILTER_DOUBLETS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_filter_doublets', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc:0.5"

    // cache false

    input:
    path archr_project
    val archr_filter_ratio

    output:
    path "proj_doublet_filtered.rds", emit: archr_project
    path "summary_filter_doublets.txt", emit: summary

    script:
    // for unknown reason, #!/usr/bin/R + direct R codes won't work
    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- filterDoublets(proj, filterRatio = $archr_filter_ratio)
    saveRDS(proj2, file = "proj_doublet_filtered.rds")

    ' > run.R

    Rscript run.R

    cp .command.log summary_filter_doublets.txt

    """
}
