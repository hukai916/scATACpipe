// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process BUILD_GENOME_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'build_genome_annotation', publish_id:'') }

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
    path bsgenome
    path gene_annotation
    path bed

    output:
    path "*.RDS", emit: genome_annotation

    // beforeScript "cp $projectDir/bin/create_geneAnnotation_genomeAnnotation.R ./"
    // also works but not best

    script:

    """
    #!/usr/bin/env Rscript

    source("$projectDir/bin/create_geneAnnotation_genomeAnnotation.R")

    bsgenome <- readRDS("$bsgenome")
    gene_annotation <- readRDS("$gene_annotation")

    if (!file.exists("$bed")) {
      stop("The supplied blacklist bed file is not valid!")
    }
    if ("$bed" == "file_token.txt") {
      bed = NULL
    } else {
      bed = "$bed"
    }

    create_ArchR_genomeannotation(
      BSgenome = bsgenome,
      geneAnnotation = gene_annotation,
      blacklist_bed = bed,
      out_dir = "./",
      $options.args)

    """
}
