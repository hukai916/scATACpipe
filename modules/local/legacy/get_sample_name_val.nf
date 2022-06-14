// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_SAMPLE_NAME_VAL {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_sample_name_val', publish_id:'') }
    container "hukai916/sinto_xenial:0.1"

    input:
    path sample

    output:
    val "${sample.baseName}", emit: sample_name_val

    script:

    """
    echo $sample > ${sample}.info.txt

    """
}
