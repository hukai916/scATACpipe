// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_GET_ANNOTATION_BIOC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_annotation_custom', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path txdb
    path org
    path bsgenome
    val archr_thread

    output:
    path "genomeAnnotation.rds", emit: genomeAnnotation
    path "geneAnnotation.rds", emit: geneAnnotation
    path "user_rlib", emit: user_rlib

    script:

    """
    #!/usr/bin/env Rscript

    library(ArchR)
    library(tools)
    library(AnnotationDbi)

    addArchRThreads(threads = $archr_thread)

    # Prepare genome/gene annotation
    dir.create("user_rlib", recursive = TRUE)
    .libPaths("user_rlib")

    ## Install BSgenome:
    if (!requireNamespace("$bsgenome", quietly = TRUE)){
      BiocManager::install("$bsgenome")
    }
    library($bsgenome)

    ## Install TxDb:
    if (!requireNamespace("$txdb", quietly = TRUE)){
      BiocManager::install("$txdb")
    }
    library($txdb)

    ## Install OrgDb:
    if (!requireNamespace("$org", quietly = TRUE)){
      BiocManager::install("$org")
    }
    library($org)

    ## Create ArchR gene/genome annotation files:
    genomeAnnotation <- createGenomeAnnotation(genome = $bsgenome)
    geneAnnotation <- createGeneAnnotation(TxDb = $txdb, OrgDb = $org)

    saveRDS(genomeAnnotation, file = "genomeAnnotation.rds")
    saveRDS(geneAnnotation, file = "geneAnnotation.rds")

    """
}
