// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_SAMPLE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'merge_sample', publish_id:'') }
    container "hukai916/miniconda3_xenial:0.1"

    input:
    val sample_name
    path sample_files

    output:
    tuple val(sample_name), path("*_merge_read1.fastq.gz"), path("*_merge_read2.fastq.gz"), path("*_merge_barcode.fastq.gz"), emit: sample_name_r1_r2_barcode
    tuple val(sample_name), path("*_merge_barcode.fastq.gz"), emit: sample_name_barcode

    script:

    """
    # merge all lanes that are from the same sample:
    cat ${sample_name}_S1_L*_R1_001.fastq.gz > ${sample_name}_merge_read1.fastq.gz
    cat ${sample_name}_S1_L*_R3_001.fastq.gz > ${sample_name}_merge_read2.fastq.gz
    cat ${sample_name}_S1_L*_R2_001.fastq.gz > ${sample_name}_merge_barcode.fastq.gz

    """

}
