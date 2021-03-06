// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_CALL_PEAKS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_call_peaks', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val groupby
    val archr_thread

    output:
    path "proj_call_peaks.rds", emit: archr_project

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes().intdiv(task.cpus)}" : ''
    // Ref: https://gitter.im/nextflow-io/nextflow?at=5a4f8f01ce68c3bc7480d7c5

    """
    echo '
    library(ArchR)
    
    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    # Add called peaks:
    pathToMacs2 <- findMacs2()
    proj2 <- addReproduciblePeakSet(
      ArchRProj = proj,
      groupBy = "$groupby",
      pathToMacs2 = pathToMacs2,
      $options.args
    )

    proj2 <- addPeakMatrix(proj2, force = TRUE)
    saveRDS(proj2, file = "proj_call_peaks.rds")

    ' > run.R

    Rscript run.R

    """
}
