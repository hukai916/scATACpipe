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
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)
    path barcode_whitelist_folder

    output:
    tuple val(sample_name), path(barcode_fastq), path("whitelist_*") emit: barcode_whitelist

    script:
    
    """
    # Use the first 10,000 barcode reads for quick test:
    # note "|| true" is to capture and skip the SIGPIPE error
    (zcat $barcode_fastq || true) | awk '{ print \$0 } NR==40000 {exit}' | gzip > subset.fastq.gz
    cat barcode_${sample_name}_*.fastq.gz > barcode_${sample_name}.fastq.gz

    get_whitelist_barcode.py $barcode_whitelist_folder subset.fastq.gz whitelist_${sample_name}_

    """
}
