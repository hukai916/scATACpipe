// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AMULET_MERGE_DOUBLETS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'amulet_merge_doublets', publish_id:'') }
    container "hukai916/amulet_xenial:0.1"

    input:
    tuple val(sample_name), path(doublets)

    output:
    tuple val(sample_name), path("merged.txt"), emit: doublets

    script:

    """
    # to do here
    touch merged.txt
    """
}
