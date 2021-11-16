// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_GET_CLUSTERING_TSV {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_clustering_tsv', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc:0.5"

    // cache false

    input:
    path archr_project
    tuple val(sample_name), path(fragment)
    val cluster

    output:
    tuple val(sample_name), path(fragment), path("*.tsv"), emit: res
    path "*.tsv", emit: tsv

    script:

    """
    echo '
    library(ArchR)
    library(stringr)

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
