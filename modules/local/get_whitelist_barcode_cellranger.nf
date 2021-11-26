// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_WHITELIST_BARCODE_CELLRANGER {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_whitelist_barcode_cellranger', publish_id:'') }
    container "hukai916/seqkit_0.16.1:0.1"

    input:
    tuple val(sample_name), path(sample_files), path(bam), path(fragment)
    path barcode_whitelist_folder

    output:
    tuple val(sample_name), path(sample_files), path(bam), path(fragment), emit: sample
    path "selected_*", emit: whitelist_barcode

    script:

    """
    # Use the first 10,000 reads for quick tests:
    barcode_files=(\$(ls -d ${sample_name}*_R2_001.fastq.gz))
    (zcat \${barcode_files[0]} || true) | awk '{ print \$0 } NR==40000 {exit}' | gzip > subset.fastq.gz

    get_whitelist_barcode.py $barcode_whitelist_folder subset.fastq.gz selected_

    """
}
