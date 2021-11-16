// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_GET_ANNOTATION_CUSTOM {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_annotation_custom', publish_id:'') }

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
    path txdb
    path org
    path bsgenome

    output:
    path "genomeAnnotation.rds", emit: genomeAnnotation
    path "geneAnnotation.rds", emit: geneAnnotation

    script:

    """
    echo '
    library(ArchR)
    library(tools)
    library(AnnotationDbi)

    # Prepare genome/gene annotation
    # TxDb must be saved/loaded with saveDb/loadDb:
    if (!(file_ext("$txdb") == "sqlite")) {
      stop("--txdb must refer to a TxDb object saved in .sqlite format!")
    } else {
      txdb <- loadDb(file = "$txdb")
    }
    if (!(file_ext("$org") == "sqlite")) {
      stop("--org must refer to a OrgDb object saved in .sqlite format!")
    } else {
      org <- loadDb(file = "$org")
    }
    if (!(tolower(file_ext("$bsgenome")) == "rds")) {
      stop("--bsgenome must refer to a BSgenome object saved in .rds/.RDS format!")
    } else {
      bsgenome <- readRDS("$bsgenome")
    }

    genomeAnnotation <- createGenomeAnnotation(genome = bsgenome)
    geneAnnotation <- createGeneAnnotation(TxDb = txdb, OrgDb = org)

    saveRDS(genomeAnnotation, file = "genomeAnnotation.rds")
    saveRDS(geneAnnotation, file = "geneAnnotation.rds")

    dir.create("user_rlib", recursive = TRUE)
    # .libPaths("user_rlib")
    # library(devtools)

    # devtools::check(package_seed[1])
    # devtools::build(package_seed[1])
    # devtools::install(package_seed[1])


    ' > run.R

    Rscript run.R

    """
}
