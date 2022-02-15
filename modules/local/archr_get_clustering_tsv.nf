// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_GET_CLUSTERING_TSV {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_clustering_tsv', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

    input:
    path archr_project
    tuple val(sample_name), path(fragment)
    val archr_thread

    output:
    tuple val(sample_name), path(fragment), path("*.tsv"), emit: res
    path "*.tsv", emit: tsv

    script:

    """
    echo '
    library(ArchR)
    library(stringr)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    index <- which(proj\$Sample == "$sample_name")
    barcodes <- str_sub(proj\$cellNames[index], nchar("$sample_name") + 2, end = -1)

    # clusters <- c("Clusters_Seurat_IterativeLSI", "Clusters_Scran_IterativeLSI", "Clusters_Seurat_Harmony", "Clusters_Scran_Harmony", "Clusters2_Seurat_IterativeLSI", "Clusters2_Scran_IterativeLSI", "Clusters2_Seurat_Harmony", "Clusters2_Scran_Harmony")
    clusters <- c("Clusters", "Clusters2")

    for (cluster in clusters) {
      tryCatch({
        eval(str2lang(paste0("clustering <- proj\$", cluster, "[index]")))
        df <- data.frame(barcodes = barcodes, clustering = clustering)
        write.table(df, file=paste0(cluster, "_", "$sample_name", ".tsv"), quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
      },
        error=function(e) {
          message(paste0("Skipping generating tsv for ", cluster, "!"))
        }
      )
    }
    ' > run.R

    Rscript run.R

    """
}
