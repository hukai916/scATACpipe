// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_WHITELIST_BARCODE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_whitelist_barcode', publish_id:'') }
    container "hukai916/seqkit_0.16.1:0.1"

    input:
    val sample_name
    path barcode_fastq
    path barcode_whitelist_folder
    // tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)
    // path barcode_whitelist_folder

    output:
    tuple val(sample_name), path("barcode_${sample_name}_combined.fastq.gz"), path("whitelist_*"), emit: sample_name_barcode_whitelist

    script:

    """
    # step1: concatenate all barcode reads or read chunks that belong to the sample_name:
    cat barcode_${sample_name}_*.fastq.gz > barcode_${sample_name}_combined.fastq.gz

    # step2: determine the whitelist barcode, use 10,000 reads for quick test:
    (zcat barcode_${sample_name}_combined.fastq.gz || true) | awk '{ print \$0 } NR==40000 {exit}' | gzip > subset.fastq.gz
    # note "|| true" is to capture and skip the SIGPIPE error

    get_whitelist_barcode.py $barcode_whitelist_folder subset.fastq.gz whitelist_${sample_name}_

    """
}
