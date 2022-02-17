// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_BATCH_CORRECTION {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_batch_correction', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

    input:
    path archr_project
    val archr_thread

    output:
    path "proj_batch_correct.rds", emit: archr_project

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    if (length(proj@sampleColData[[1]]) < 2) {
      stop("Only 1 sample detected for Harmony!")
    }

    proj2 <- addHarmony(
      ArchRProj = proj,
      reducedDims = "IterativeLSI",
      name = "Harmony",
      groupBy = "Sample",
      $options.args
    )

    saveRDS(proj2, file = "proj_batch_correct.rds")
    ' > run.R

    Rscript run.R

    """
}
