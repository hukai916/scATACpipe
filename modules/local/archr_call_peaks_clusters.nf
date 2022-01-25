// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_CALL_PEAKS_CLUSTERS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_call_peaks_clusters', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    path user_rlib
    val archr_thread

    output:
    path "proj_call_peaks.rds", emit: archr_project

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes().intdiv(task.cpus)}" : ''
    // Ref: https://gitter.im/nextflow-io/nextflow?at=5a4f8f01ce68c3bc7480d7c5

    """
    echo '
    library(ArchR)
    .libPaths("user_rlib") # for user installed packages

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    # Add called peaks:
    pathToMacs2 <- findMacs2()
    clusters <- c("Clusters_Seurat_IterativeLSI", "Clusters_Scran_IterativeLSI", "Clusters_Seurat_Harmony", "Clusters_Scran_Harmony", "Clusters2_Seurat_IterativeLSI", "Clusters2_Scran_IterativeLSI", "Clusters2_Seurat_Harmony", "Clusters2_Scran_Harmony")
    for (cluster in clusters) {
      tryCatch({
        proj <- addReproduciblePeakSet(
          ArchRProj = proj,
          groupBy = cluster,
          pathToMacs2 = pathToMacs2,
          $options.args
        )
        proj <- addPeakMatrix(proj, force = TRUE)
      },
        error=function(e) {
          message(paste0("Skipping calling peaks for ", cluster, "!"))
        }
      )
    }

    saveRDS(proj, file = "proj_call_peaks.rds")

    ' > run.R

    Rscript run.R

    """
}
