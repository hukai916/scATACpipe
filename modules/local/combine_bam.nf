// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COMBINE_BAM {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'combine_bam', publish_id:'') }
    container "hukai916/pysam_xenial:0.1"

    input:
    val sample_name
    path bam

    output:
    val sample_name, emit: sample_name
    path "*.combined.sorted.bam", emit: bam

    script:

    """
    # first combine all bam files that belong to the same library (sample_name):
    samtools merge ${sample_name}.combined.bam ${sample_name}.*.bam

    # then, sort output bam:
    samtools sort ${sample_name}.combined.bam -o ${sample_name}.combined.sorted.bam

    """

}
