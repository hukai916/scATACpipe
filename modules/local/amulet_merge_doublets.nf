// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AMULET_MERGE_DOUBLETS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'amulet_merge_doublets', publish_id:'') }
    container "hukai916/miniconda3_v24.7_amd64_bio:0.1"

    input:
    path doublets

    output:
    path "amulet_doublets.txt", emit: amulet_doublets

    script:

    """
    cat cells_filter_*.txt > amulet_doublets.txt

    """
}
