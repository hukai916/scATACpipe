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


    script:

    """
    # filter fragment file
    # filter_fragment.py $fragment $valid_barcode | gzip > valid_barcode_filtered_fragment.gz

    # filter bam file

    """
}
