// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FRAG_TO_FREQ {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'frag_to_freq', publish_id:'') }
    container "hukai916/pysam_xenial:0.1"

    input:
    tuple val(sample_name), path(fragments)

    output:
    tuple val(sample_name), path(fragments), path("freq_${sample_name}.sorted.tsv.gz"), emit: sample_name_frag_freq

    script:

    """
    echo $options.args > test.txt
    frag_to_freq.py $fragments freq_${sample_name}.sorted.tsv
    bgzip freq_${sample_name}.sorted.tsv

    """

}
