// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process BUILD_GENE_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'build_gene_annotation', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_util:0.1"

    // cache false

    input:
    path txdb
    path gtf
    val species_latin_name

    output:
    path "*", emit: gene_annotation

    // beforeScript "cp $projectDir/bin/create_geneAnnotation_genomeAnnotation.R ./"
    // also works but not best

    script:

    """
    #!/usr/bin/env Rscript

    source("$projectDir/bin/create_geneAnnotation_genomeAnnotation.R")
    txdb <- loadDb("$txdb")

    # first, get gene symbol SimpleList
    if (is.null($species_latin_name)) {
      id2symbol <- get_geneID_symbol(gtf = "$gtf", species_latin_name = NULL)
    } else {
      id2symbol <- get_geneID_symbol(gtf = "$gtf", species_latin_name = "$species_latin_name")
    }

    create_ArchR_geneannotation_WO_OrgDb(TxDb = txdb,
                                  geneID2Symbol = id2symbol,
                                  out_dir = "./",
                                  $options.args
                                  )

    """
}
