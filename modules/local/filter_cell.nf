// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTER_CELL {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'filter_cell', publish_id:'') }

    container "hukai916/sinto_xenial:0.1"

    input:
    tuple val(sample_name), path(sample_files), path(bam), path(fragment)
    path valid_barcode

    output:
    val sample_name, emit: sample_name
    path "*_valid_barcode_filtered_fragment.tsv.gz", emit: filtered_fragment
    tuple val(sample_name), path("*_valid_barcode_filtered_fragment.tsv.gz"), emit: ch_filtered_fragment
    path "*_valid_barcode_filtered_fragment.bam", emit: filtered_bam

    script:

    """
    # filter fragment file
    filter_fragment.py $fragment $valid_barcode | gzip > ${sample_name}_valid_barcode_filtered_fragment.tsv.gz

    # filter bam file
    samtools index $bam
    filter_bam.py $bam $valid_barcode ${sample_name}_valid_barcode_filtered_fragment.bam

    """
}
