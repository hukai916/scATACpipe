// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_FRAGMENTS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'fragments', publish_id:'') }
    container "hukai916/sinto_xenial:0.2"

    input:
    tuple val(sample_name), path(bam)

    output:
    val sample_name, emit: sample_name
    path "*.sorted.tsv.gz", emit: fragments
    tuple val(sample_name), path("*.sorted.tsv.gz"), emit: ch_fragment

    script:

    """
    # first index the bam file
    samtools index $options.args $bam

    # then, generate the fragments file
    sinto fragments $options.args --nproc $task.cpus --bam $bam -f fragments.tsv --barcode_regex "[^:]*"
    # sort the fragment, note that bgzip is a must
    sort -k 1,1 -k2,2n fragments.tsv | bgzip > ${sample_name}.sorted.tsv.gz

    # clean temp file:
    rm fragments.tsv

    """
}
