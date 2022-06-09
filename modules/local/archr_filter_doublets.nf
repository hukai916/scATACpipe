// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_FILTER_DOUBLETS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_filter_doublets', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path archr_project
    val archr_thread

    output:
    path "proj_doublet_filtered.rds", emit: archr_project
    path "summary_filter_doublets.txt", emit: summary

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- filterDoublets(proj, filterRatio = $options.archr_filter_doublets_ratio)
    saveRDS(proj2, file = "proj_doublet_filtered.rds")

    ' > run.R

    Rscript run.R

    cp .command.log summary_filter_doublets.txt

    """
}
