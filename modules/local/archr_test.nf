// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_TEST {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_test', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project

    output:
    path "file.txt", emit: test_file

    script:

    """
    echo $archr_project > file.txt
    """
}
