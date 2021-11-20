// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_PSEUDO_BULK {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_pseudo_bulk', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val groupby

    output:
    path "proj_pseudo_bulk.rds", emit: archr_project

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    # Add pseudo-bulk:
    proj2 <- addGroupCoverages(ArchRProj = proj, groupBy = "$groupby")
    saveRDS(proj2, file = "proj_pseudo_bulk.rds")

    ' > run.R

    Rscript run.R

    """
}
