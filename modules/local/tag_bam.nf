// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TAG_BAM {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'tag_bam', publish_id:'') }

    container "hukai916/sinto_xenial:0.2"

    input:
    tuple val(sample_name), val(chunk_name), path(tagfile), path(bam)
    // val sample_name
    // path tagfile
    // path bam

    output:
    val sample_name, emit: sample_name
    path "*.tag.bam", emit: bam

    script:

    """
    samtools index $bam
    tag_bam.py $bam $tagfile ${chunk_name}.tag.bam $task.cpus

    """
}
