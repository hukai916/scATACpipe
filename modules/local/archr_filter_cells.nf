// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_FILTER_CELLS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_filter_cells', publish_id:'') }
    container "hukai916/r_env_v4.4.1_amd64:0.1"

    input:
    path archr_project
    val archr_thread

    output:
    path "proj_cell_filtered.rds", emit: archr_project

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    filters <- trimws(strsplit("$options.args", split = ",")[[1]], "both")
    idxPass <- c()
    for (filter in filters) {
      command <- paste0("idxPass_tem <- which(proj\$", filter, ")")
      print(command)
      eval(str2lang(command))
      if (length(idxPass) == 0) {
        idxPass <- idxPass_tem
      } else {
        idxPass <- intersect(idxPass, idxPass_tem)
      }
    }

    cellsPass <- proj\$cellNames[idxPass]
    proj2     <- proj[cellsPass,]

    saveRDS(proj2, file = "proj_cell_filtered.rds")

    ' > run.R

    Rscript run.R

    """
}
