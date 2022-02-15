// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_CALL_PEAKS_CLUSTERS2 {
  // Note that each ArchRProj can only have one set of called peaks
  // Duplicate the ArchRProj at pseudo_bulk step already.

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_call_peaks_clusters2', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

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
    library(parallel) # quick fix for Haibo ArchR
    .libPaths("user_rlib") # for user installed packages

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    pathToMacs2 <- findMacs2()
    eval(str2lang(paste0("library(", getGenome(proj), ")")))
    genomeSizeCmd <- paste0("sum(", getGenome(proj), "@seqinfo@seqlengths)")
    genomeSize <- eval(str2lang(genomeSizeCmd))

    proj <- addReproduciblePeakSet(
      ArchRProj = proj,
      groupBy = "Clusters2",
      pathToMacs2 = pathToMacs2,
      genomeSize = genomeSize
      $options.args
    )
    proj <- addPeakMatrix(proj, force = TRUE)

    saveRDS(proj, file = "proj_call_peaks.rds")

    ' > run.R

    Rscript run.R

    """
}
