// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_PSEUDO_BULK_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_pseudo_bulk_clusters', publish_id:'') }

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

    output:
    path "archr_project.rds", emit: archr_project
    path "save_archr_project", emit: archr_dir

    script:

    """
    echo '
    library(ArchR)

    proj <- readRDS("$archr_project", refhook = NULL)
    proj2 <- saveArchRProject(ArchRProj = proj, outputDirectory = "save_archr_project", load = TRUE)

    # Add pseudo-bulk:
    # addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters", force = TRUE)
    # Tested, above will add GroupCoverage, but wont update proj2 automatically, must use below
    proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters", force = TRUE)
    saveRDS(proj2, file = "archr_project.rds")

    ' > run.R

    Rscript run.R

    """
}
