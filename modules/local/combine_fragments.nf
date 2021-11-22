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
    path "combined.fragments.sort.bed.gz", emit: fragments
    tuple val(sample_name), path("combined.fragments.sort.bed.gz"), emit: ch_fragment

    script:

    """
    # first combine all fragments that belong to the same library (sample_name)


    # first index the bam file
    samtools index $options.args $bam

    # then, generate the fragments file
    sinto fragments $options.args --nproc $task.cpus --bam $bam -f fragments.bed --barcode_regex "[^:]*"
    # sort and bzip the fragment file
    sort -k 1,1 -k2,2n fragments.bed > fragments.sort.bed
    bgzip fragments.sort.bed

    """
}
