// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_CALL_PEAKS_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_call_peaks_clusters', publish_id:'') }

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
    path "proj_call_peaks.rds", emit: archr_project

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes().intdiv(task.cpus)}" : ''
    // Ref: https://gitter.im/nextflow-io/nextflow?at=5a4f8f01ce68c3bc7480d7c5

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    # Add called peaks:
    pathToMacs2 <- findMacs2()
    proj2 <- addReproduciblePeakSet(
      ArchRProj = proj,
      groupBy = "Clusters",
      pathToMacs2 = pathToMacs2,
      $options.args
    )

    proj2 <- addPeakMatrix(proj2, force = TRUE)
    saveRDS(proj2, file = "proj_call_peaks.rds")

    ' > run.R

    Rscript run.R

    """
}
