// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUTADAPT {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'cutadapt', publish_id:'') }
    container "hukai916/cutadapt_xenial:0.1"

    input:
    val sample_name
    path read1_fastq
    path read2_fastq
    val read1_adapter
    val read2_adapter

    output:
    val sample_name, emit: sample_name
    path "R1/trimmed*", emit: trimed_read1_fastq
    path "R2/trimmed*", emit: trimed_read2_fastq
    path "log_cutadapt_*.txt", emit: log

    script:
    read1_name = read1_fastq.getName()
    read2_name = read2_fastq.getName()

    """
    mkdir R1
    cutadapt $options.args -a $read1_adapter -o R1/trimmed_$read1_name $read1_fastq
    mkdir R2
    cutadapt $options.args -a $read2_adapter -o R2/trimmed_$read2_name $read2_fastq
    cp .command.log log_cutadapt_${sample_name}.txt

    """
}
