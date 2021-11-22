// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MATCH_READS_TRIMMED {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'match_reads_trimmed', publish_id:'') }
    container "hukai916/fastq-pair:0.1"

    input:
    val sample_name
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name
    path "match_trim/first_read_in_pair.fq.paired.fq.gz", emit: read1_fastq
    path "match_trim/second_read_in_pair.fq.paired.fq.gz", emit: read2_fastq

    script:

    """
    mkdir match_trim
    cp $read1_fastq match_trim/
    cp $read2_fastq match_trim/
    cd match_trim
    mv $read1_fastq first_read_in_pair.fq.gz
    mv $read2_fastq second_read_in_pair.fq.gz
    gzip -d first_read_in_pair.fq.gz
    gzip -d second_read_in_pair.fq.gz

    fastq_pair $options.args first_read_in_pair.fq second_read_in_pair.fq
    rm first_read_in_pair.fq second_read_in_pair.fq
    gzip first_read_in_pair.fq.paired.fq
    gzip second_read_in_pair.fq.paired.fq

    """
}
