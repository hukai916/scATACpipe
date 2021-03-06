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
    path barcode_chunk
    path barcode_whitelist_folder

    output:
    // tuple val(sample_name), path("barcode_*.fastq.gz"), path("R1_*.fastq.gz"), emit: sample
    val sample_name, emit: sample_name
    path "whitelist_*", emit: whitelist_barcode
    path "barcode_*.fastq.gz", emit: barcode_fastq
    // path "R1_*.fastq.gz", emit: read1_fastq

    script:

    """
    # Use the first 10,000 reads from one chunk for quick tests:
    barcode_files=(\$(ls -d barcode_${sample_name}_*.fastq.gz))

    # note "|| true" is to capture and skip the SIGPIPE error
    (zcat \${barcode_files[0]} || true) | awk '{ print \$0 } NR==40000 {exit}' | gzip > subset.fastq.gz
    cat barcode_${sample_name}_*.fastq.gz > barcode_${sample_name}.fastq.gz
    #cat R1_${sample_name}_*.fastq.gz > R1_${sample_name}.fastq.gz

    get_whitelist_barcode.py $barcode_whitelist_folder subset.fastq.gz whitelist_${sample_name}_

    """
}
