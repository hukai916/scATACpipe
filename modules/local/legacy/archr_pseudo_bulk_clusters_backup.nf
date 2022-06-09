// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_PSEUDO_BULK_CLUSTERS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_pseudo_bulk_clusters', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

    input:
    path archr_project
    path user_rlib
    val archr_thread

    output:
    path "archr_project_pseudobulk.rds", emit: archr_project
    path user_rlib, emit: user_rlib

    script:

    """
    echo '
    library(ArchR)
    .libPaths("user_rlib") # for user installed packages

    addArchRThreads(threads = $archr_thread)

    proj  <- readRDS("$archr_project", refhook = NULL)

    clusters <- c("Clusters_Seurat_IterativeLSI", "Clusters_Scran_IterativeLSI", "Clusters_Seurat_Harmony", "Clusters_Scran_Harmony", "Clusters2_Seurat_IterativeLSI", "Clusters2_Scran_IterativeLSI", "Clusters2_Seurat_Harmony", "Clusters2_Scran_Harmony")
    for (cluster in clusters) {
      tryCatch({
        proj <- addGroupCoverages(ArchRProj = proj, groupBy = cluster, force = TRUE)
      },
        error=function(e) {
          message(paste0("Skipping adding pseudo-bulk for ", cluster, "!"))
        }
      )
    }

    saveRDS(proj, file = "archr_project_pseudobulk.rds")

    ' > run.R

    Rscript run.R

    """
}
