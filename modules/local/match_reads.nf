// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MATCH_READS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'match_reads', publish_id:'') }
    container "hukai916/seqkit_0.16.1:0.1"

    input:
    val sample_name
    path corrected_barcode_fastq
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name
    path corrected_barcode_fastq, emit: barcode_fastq
    path "R1/*.fastq.gz", emit: read1_fastq
    path "R2/*.fastq.gz", emit: read2_fastq
    path "R1_barcode/*.fastq.gz", emit: barcode1_fastq
    path "R2_barcode/*.fastq.gz", emit: barcode2_fastq

    script:

    """
    seqkit pair $options.args -1 $corrected_barcode_fastq -2 $read1_fastq -O R1
    seqkit pair $options.args -1 $corrected_barcode_fastq -2 $read2_fastq -O R2

    mkdir R1_barcode R2_barcode
    mv R1/$corrected_barcode_fastq R1_barcode/$corrected_barcode_fastq
    mv R2/$corrected_barcode_fastq R2_barcode/$corrected_barcode_fastq

    # Note, by using this method here, R1 paired up with R1 barcode and R2 paired up with R2 barcode, but R1 may not pair up with R1.
    # It will not hurt since the match_reads_trimmed.nf pairs up R1 and R2 after trimming.


    """
}
