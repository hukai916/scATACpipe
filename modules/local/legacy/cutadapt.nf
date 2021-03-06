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
    tuple val(sample_name), path(read1_fastq), path(read2_fastq)
    val read1_adapter
    val read2_adapter

    output:
    tuple val(sample_name), path("R1/*trimmed.fastq.gz"), path("R2/*trimmed.fastq.gz"), emit: reads_0
    path "log_cutadapt_*.txt", emit: log

    script:
    read1_trimmed_fastq = read1_fastq.split("\\.")[0..-3].join(".") + ".trimmed.fastq.gz"
    read2_trimmed_fastq = read2_fastq.split("\\.")[0..-3].join(".") + ".trimmed.fastq.gz"


    """
    sample_name=$read1_fastq
    outname1="\${sample_name%%.*}"
    sample_name=$read2_fastq
    outname2="\${sample_name%%.*}"
    outname="\${outname1#R1_}"

    mkdir R1 R2
    cutadapt $options.args -a $read1_adapter -A $read2_adapter -o R1/\${outname1}.trimmed.fastq.gz -p R2/\${outname2}.trimmed.fastq.gz $read1_fastq $read2_fastq
    cp .command.log log_cutadapt_\${outname}.txt

    """
}
