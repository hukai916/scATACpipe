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
    tuple val(sample_name), path(read1_fastq), path(read2_fastq)

    output:
    tuple val(sample_name), path("match_trim/*first_read_in_pair.fq.paired.fq.gz"), path("match_trim/*second_read_in_pair.fq.paired.fq.gz"), emit: reads_0

    script:

    """
    mkdir match_trim
    cp $read1_fastq match_trim/ # if using cp -P, gzip is problematic then.
    cp $read2_fastq match_trim/
    cd match_trim
    mv $read1_fastq ${sample_name}_first_read_in_pair.fq.gz
    mv $read2_fastq ${sample_name}_second_read_in_pair.fq.gz
    gzip -d ${sample_name}_first_read_in_pair.fq.gz
    gzip -d ${sample_name}_second_read_in_pair.fq.gz

    fastq_pair $options.args ${sample_name}_first_read_in_pair.fq ${sample_name}_second_read_in_pair.fq
    rm ${sample_name}_first_read_in_pair.fq ${sample_name}_second_read_in_pair.fq
    gzip ${sample_name}_first_read_in_pair.fq.paired.fq
    gzip ${sample_name}_second_read_in_pair.fq.paired.fq

    """
}
