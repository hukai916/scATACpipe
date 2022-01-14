// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AMULET_FILTER_DOUBLETS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'amulet_filter_doublets', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    path cells_filter

    output:
    path "proj_doublet_filtered.rds", emit: archr_project
    // path "summary_filter_doublets.txt", emit: summary
    tuple path(archr_project), path(cells_filter), emit: test_input

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    cellsFilter <- scan("$cells_filter", what = "character")
    proj@cellColData <- proj@cellColData[rownames(proj@cellColData) %ni% cellsFilter,,drop=FALSE]

    saveRDS(proj, file = "proj_doublet_filtered.rds")

    ' > run.R

    Rscript run.R

    cp .command.log summary_amulet_filter_doublets.txt
    """
}
