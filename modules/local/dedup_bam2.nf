// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DEDUP_BAM2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'dedup_bam', publish_id:'') }
    container "hukai916/sinto_kai:0.3"

    input:
    val sample_name
    path bam
    val barcode_tag

    output:
    tuple val(sample_name), path("*.dedup.bam"), emit: sample_name_bam
    path "*.dedup.bam", emit: bam
    path "summary_rm_*", emit: dedup_summary

    script:

    """
    # Deduplicate bam file with remove_duplicate.py:
    remove_duplicate2.py --inbam $bam --barcode_tag $barcode_tag --extend_softclip 0 --outdir ./ --outbam ${sample_name}.dedup.bam --nproc $task.cpus

    """

}
