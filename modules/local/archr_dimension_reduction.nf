// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_DIMENSION_REDUCTION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_dimension_reduction', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project

    output:
    path "proj_lsi.rds", emit: archr_project

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- addIterativeLSI(
      ArchRProj = proj,
      useMatrix = "TileMatrix",
      name = "IterativeLSI",
      $options.args
    )

    saveRDS(proj2, file = "proj_lsi.rds")
    ' > run.R

    Rscript run.R

    """
}
