// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DEDUP_BAM {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'dedup_bam', publish_id:'') }
    container "hukai916/miniconda3_picard:0.1"

    input:
    val sample_name
    path bam

    output:
    val sample_name, emit: sample_name
    path "*.dedup.bam", emit: bam

    script:

    """
    # Deduplicate bam file with Picard::markDuplicates:
    picard MarkDuplicates REMOVE_DUPLICATES=true I=$bam O=${sample_name}.dedup.bam M=matrix_${sample_name}.txt

    """

}
