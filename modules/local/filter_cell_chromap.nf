// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTER_CELL_CHROMAP {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'filter_cell_chromap', publish_id:'') }

    container "hukai916/miniconda3_v24.7_amd64_bio:0.1"

    input:
    tuple val(sample_name), path(fragment), path(filtered_barcode)

    output:
    val sample_name, emit: sample_name
    path "*_valid_barcode_filtered_fragment.tsv.gz", emit: filtered_fragment
    tuple val(sample_name), path("*_valid_barcode_filtered_fragment.tsv.gz"), emit: sample_name_filtered_fragment

    script:

    """
    # filter fragment file
    filter_fragment.py $fragment $filtered_barcode | bgzip > ${sample_name}_valid_barcode_filtered_fragment.tsv.gz

    """
}
