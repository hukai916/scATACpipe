// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTER_CELL {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'filter_cell', publish_id:'') }

    container "hukai916/sinto_xenial:0.2"

    input:
    tuple val(sample_name), path(bam), path(fragment), path(filtered_barcode)
    val(barcode_tag)

    output:
    val sample_name, emit: sample_name
    path "*_valid_barcode_filtered_fragment.tsv.gz", emit: filtered_fragment
    tuple val(sample_name), path("*_valid_barcode_filtered_fragment.tsv.gz"), emit: sample_name_filtered_fragment
    path "*_valid_barcode_filtered_fragment.bam", emit: filtered_bam

    script:

    """
    # filter fragment file
    filter_fragment.py $fragment $filtered_barcode | bgzip > ${sample_name}_valid_barcode_filtered_fragment.tsv.gz

    # filter bam file
    samtools index $bam
    filter_bam.py $bam $filtered_barcode ${sample_name}_valid_barcode_filtered_fragment.bam $barcode_tag

    """
}
