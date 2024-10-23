// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLIT_FASTQ {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'split_fastq', publish_id:'') }
    container "hukai916/miniconda3_v24.7_amd64_bio:0.1"

    input:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)
    val sample_count

    output:
    val sample_name
    path "R1_*.fastq.gz", emit: read1_fastq
    path "R2_*.fastq.gz", emit: read2_fastq
    path "barcode_*.fastq.gz", emit: barcode_fastq

    script:

    """
    # default to 20M reads (80000000) per chunk
    pigz -dc -p $task.cpus $read1_fastq | split --lines=80000000 --filter='pigz -p $task.cpus > \${FILE}.fastq.gz' - R1_${sample_name}_${sample_count}_ &
    pigz -dc -p $task.cpus $read2_fastq | split --lines=80000000 --filter='pigz -p $task.cpus > \${FILE}.fastq.gz' - R2_${sample_name}_${sample_count}_ &
    pigz -dc -p $task.cpus $barcode_fastq | split --lines=80000000 --filter='pigz -p $task.cpus > \${FILE}.fastq.gz' - barcode_${sample_name}_${sample_count}_ &

    """
}
