// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADD_BARCODE_TO_READS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'add_barcode_to_reads', publish_id:'') }
    container "hukai916/sinto_xenial:0.2"

    input:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)
    val sample_count

    output:
    tuple val(sample_name), path("R1_${sample_name}_${sample_count}.barcoded.fastq.gz"), path("R2_${sample_name}_${sample_count}.barcoded.fastq.gz"), emit: reads_0
    tuple val(sample_name), path("R1_${sample_name}_${sample_count}.barcoded.fastq.gz"), path("R2_${sample_name}_${sample_count}.barcoded.fastq.gz"), path("barcode_${sample_name}_${sample_count}.fastq.gz"), emit: reads
    // Below are for GET_WHITELIST_BARCODE
    val sample_name, emit: sample_name
    path "barcode_${sample_name}_${sample_count}.fastq.gz", emit: barcode_fastq

    script:
    read1_barcoded_fastq = read1_fastq.name.split("\\.")[0..-3].join(".") + ".barcoded.fastq.gz"
    read2_barcoded_fastq = read2_fastq.name.split("\\.")[0..-3].join(".") + ".barcoded.fastq.gz"

    """
    # use the first read length from fastq file to determine the length since -b is required by sinto.
    barcode_length=\$((zcat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    # sinto barcode $options.args --barcode_fastq $barcode_fastq --read1 $read1_fastq --read2 $read2_fastq -b \$barcode_length
    addbarcodes_parallel.py $barcode_fastq \$barcode_length $read1_fastq $read2_fastq $task.cpus

    # remove sequence description in + line, otherwise, cutadapt may complain if it does not match with line 1:
    zcat $read1_barcoded_fastq | awk '{if (\$1 ~/^\+/) {print "+"} else {print \$0}}' | gzip > tem_$read1_barcoded_fastq
    zcat $read2_barcoded_fastq | awk '{if (\$1 ~/^\+/) {print "+"} else {print \$0}}' | gzip > tem_$read2_barcoded_fastq
    mv tem_$read1_barcoded_fastq $read1_barcoded_fastq
    mv tem_$read2_barcoded_fastq $read2_barcoded_fastq

    # rename the files:
    mv $read1_barcoded_fastq R1_${sample_name}_${sample_count}.barcoded.fastq.gz
    mv $read2_barcoded_fastq R2_${sample_name}_${sample_count}.barcoded.fastq.gz
    mv $barcode_fastq barcode_${sample_name}_${sample_count}.fastq.gz # rename input is fine to nextflow

    """
}
