// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_WHITELIST_CHROMAP {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_whitelist_chromap', publish_id:'') }
    container "hukai916/seqkit_0.16.1:0.1"

    input:
    tuple val(sample_name), path(barcode_fastq)
    path barcode_whitelist_folder

    output:
    tuple val(sample_name), path("whitelist_*"), emit: sample_name_whitelist

    script:

    """
    # determine the whitelist barcode, use 10,000 reads for quick test:
    (zcat $barcode_fastq || true) | awk '{ print \$0 } NR==40000 {exit}' | gzip > subset.fastq.gz
    # note "|| true" is to capture and skip the SIGPIPE error

    get_whitelist_barcode.py $barcode_whitelist_folder subset.fastq.gz whitelist_${sample_name}_

    """
}
