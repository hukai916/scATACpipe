// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MATCH_SAMPLE_NAME {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'match_sample_name', publish_id:'') }
    container "hukai916/sinto_xenial:0.2"

    input:
    tuple val(sample_name), val(path_fastq_1), val(path_fastq_2), val(path_barcode)
    val sample_count

    output:
    val sample_name, emit: sample_name
    path "name_matched_fastq/*.fastq.gz", emit: sample_files

    // path "*_R1_001.fastq.gz", emit: read1_fastq
    // path "*_R2_001.fastq.gz", emit: barcode_fastq
    // path "*_R3_001.fastq.gz", emit: read2_fastq

    script:

    """
    # rename sample names in case of unexpected inconsistency:
    # note that the lane number may not match with original, but it will not hurt anything.

    mkdir name_matched_fastq
    cp -P $path_fastq_1 name_matched_fastq/${sample_name}_S1_L${sample_count}_R1_001.fastq.gz
    cp -P $path_fastq_2 name_matched_fastq/${sample_name}_S1_L${sample_count}_R3_001.fastq.gz
    cp -P $path_barcode name_matched_fastq/${sample_name}_S1_L${sample_count}_R2_001.fastq.gz

    """
}
