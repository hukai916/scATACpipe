// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUILD_GENE_ANNOTATION {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'build_gene_annotation', publish_id:'') }
    container "hukai916/r_env_v4.4.1_amd64:0.1"

    input:
    path txdb
    path gtf
    val species_latin_name

    output:
    path "*", emit: gene_annotation

    script:

    """
    #!/usr/bin/env Rscript

    source("$projectDir/bin/create_geneAnnotation_genomeAnnotation.R")
    txdb <- loadDb("$txdb")

    # first, get gene symbol SimpleList
    id2symbol <- get_geneID_symbol(gtf = "$gtf", species_latin_name = "$species_latin_name")

    # use Strategy 2: only keep unique gene symbols, details see id2symbol in create_ArchR_geneannotation_WO_OrgDb
    create_ArchR_geneannotation_WO_OrgDb(TxDb = txdb,
                                  geneID2Symbol = id2symbol,
                                  out_dir = "./",
                                  $options.args
                                  )

    """
}
