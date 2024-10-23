// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUILD_BSGENOME {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'build_bsgenome', publish_id:'') }
    container "hukai916/r_env_v4.4.1_amd64:0.1"

    input:
    path genome_fasta

    output:
    path "custom_BSgenome.rds", emit: bsgenome
    path "user_rlib", emit: user_rlib

    script:

    """
    build_bsgenome.R --genome_fasta $genome_fasta

    """
}
