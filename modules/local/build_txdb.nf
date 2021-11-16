// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process BUILD_TXDB {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'build_txdb', publish_id:'') }

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
    path gtf

    output:
    path "*.TxDb.sqlite", emit: txdb

    script:

    """
    # First add "chr" to GTF since "chr" is added to BSgenome but not for GTF.
    # This step is no longer needed given the PREP_GTF module, but it won't hurt.
    if [[ $gtf == *.gz ]]; then
      zcat $gtf | awk 'BEGIN{FS=OFS="\\t"}/^#/{print; next} NF==9{\$1 = gensub(/(chr)?(.+)/, "chr\\\\2", "g", \$1); print}' | gzip > chrPrefixed.gtf.gz
    else
      cat $gtf | awk 'BEGIN{FS=OFS="\\t"}/^#/{print; next} NF==9{\$1 = gensub(/(chr)?(.+)/, "chr\\\\2", "g", \$1); print}' | gzip > chrPrefixed.gtf.gz
    fi

    build_txdb.R --bsgenome $bsgenome --gtf chrPrefixed.gtf.gz

    """
}
