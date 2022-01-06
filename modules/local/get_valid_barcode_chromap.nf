// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_VALID_BARCODE_CHROMAP {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_valid_barcode_chromap', publish_id:'') }
    container "hukai916/r_util:0.2"

    input:
    tuple val(sample_name), path(fragments)

    output:
    tuple val(sample_name), path("*_valid_barcodes.txt"), emit: sample_name_valid_barcodes

    script:

    """
    # output valida barcodes:
    get_valid_barcode_inflection.R --freq $fragments --outfile ${sample_name}_valid_barcodes.txt

    """
}
