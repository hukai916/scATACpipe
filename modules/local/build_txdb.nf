// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUILD_TXDB {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'build_txdb', publish_id:'') }
    container "hukai916/r_utils:0.1"

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
