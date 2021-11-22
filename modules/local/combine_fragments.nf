// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COMBINE_FRAGMENTS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'combine_fragments', publish_id:'') }
    container "hukai916/sinto_xenial:0.2"

    input:
    val sample_name
    path fragments

    output:
    val sample_name, emit: sample_name
    path ".combined.fragments.sort.tsv.gz", emit: fragments
    tuple val(sample_name), path(".combined.fragments.sort.tsv.gz"), emit: ch_fragment

    script:

    """
    # first combine all fragments that belong to the same library (sample_name):
    cat ${sample_name}* > ${sample_name}.combined.fragments.tsv

    # then sort and bgzip:
    sort -k 1,1 -k2,2n ${sample_name}.combined.fragments.tsv > ${sample_name}.combined.fragments.sort.tsv
    bgzip ${sample_name}.combined.fragments.sort.tsv

    """
}
