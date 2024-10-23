// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMPLESHEET_CHECK_FASTQ {
    tag "$samplesheet"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    container "hukai916/miniconda3_v24.7_amd64_bio:0.1"

    input:
    path samplesheet

    output:
    path '*.csv', emit: csv
    path '*.txt', emit: count

    script:

    """
    check_samplesheet_fastq.py $samplesheet samplesheet.valid.csv
    awk '{ if (NR == 1) print "sample_count"; else print NR-1 }' samplesheet.valid.csv > sample_count.txt

    """
}
