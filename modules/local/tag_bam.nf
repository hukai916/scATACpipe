// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TAG_BAM {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'tag_bam', publish_id:'') }

    container "hukai916/sinto_xenial:0.1"

    input:
    val sample_name
    path tagfile
    path dedup_bam

    output:
    val sample_name, emit: sample_name
    path "*.tag.bam", emit: bam

    script:

    """
    samtools index ${sample_name}.dedup.bam
    tag_bam.py ${sample_name}.dedup.bam $tagfile ${sample_name}.tag.bam

    """
}
