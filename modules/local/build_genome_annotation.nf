// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUILD_GENOME_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'build_genome_annotation', publish_id:'') }
    container "hukai916/r_util:0.1"

    input:
    path bsgenome
    path gene_annotation
    path bed

    output:
    path "*.RDS", emit: genome_annotation

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
