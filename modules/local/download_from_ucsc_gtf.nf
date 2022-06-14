// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DOWNLOAD_FROM_UCSC_GTF {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_from_ucsc_gtf', publish_id:'') }
    container "hukai916/miniconda3_xenial:0.1"

    input:
    val genome_name
    path dict_json

    output:
    path "*.gtf.gz", emit: gtf

    script:

    """
    gtf_link=\$(get_download_url.py $dict_json $genome_name gtf)
    wget --user-agent=" Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36" \$gtf_link -o logfile.gtf.txt

    """
}
