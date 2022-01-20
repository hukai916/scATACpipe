// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_GET_CLUSTERING_TSV {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_clustering_tsv', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    tuple val(sample_name), path(fragment)
    val cluster
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
    clustering <- proj\$Clusters[index]

    df <- data.frame(barcodes = barcodes, clustering = clustering)
    write.table(df, file=paste0("Clusters_", "$sample_name", ".tsv"), quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

    if ("$cluster" == "Clusters2") {
      clustering <- proj\$Clusters2[index]

      df <- data.frame(barcodes = barcodes, clustering = clustering)
      write.table(df, file=paste0("Clusters2_", "$sample_name", ".tsv"), quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
    }

    ' > run.R

    Rscript run.R

    """
}
