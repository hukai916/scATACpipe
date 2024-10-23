// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_PSEUDO_BULK_CLUSTERS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_pseudo_bulk_clusters', publish_id:'') }
    container "hukai916/r_env_v4.4.1_amd64:0.1"

    input:
    path archr_project
    path user_rlib
    val archr_thread

    output:
    path "archr_project_pseudobulk.rds", emit: archr_project
    path "archr_project_pseudobulk", emit: archr_dir
    path user_rlib, emit: user_rlib

    script:

    """
    echo '
    library(ArchR)
    .libPaths("user_rlib") # for user installed packages

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)
    proj2 <- saveArchRProject(ArchRProj = proj, outputDirectory = "archr_project_pseudobulk", load = TRUE)
    # Note that saveArchRProject is a must here since each ArchRProj can only have one set of PeakSet, duplicate the proj for cluster2

    proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters", force = TRUE)
    proj2 <- addImputeWeights(proj2)

    saveRDS(proj2, file = "archr_project_pseudobulk.rds")

    ' > run.R

    Rscript run.R

    """
}
