// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MATCH_READS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'match_reads', publish_id:'') }
    container "hukai916/fastq-pair:0.1"

    input:
    val sample_name
    path corrected_barcode_fastq
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name

    path "match_pair_first_read/first_read_in_pair.fq.paired.fq", emit: read1_fastq
    path "match_pair_second_read/second_read_in_pair.fq.paired.fq", emit: read2_fastq
    path "match_pair_first_read/first_read_matched_corrected_barcode.fq.gz", emit: barcode1_fastq
    path "match_pair_second_read/second_read_matched_corrected_barcode.fq.gz", emit: barcode2_fastq
    // Be carefule of the duplicated staged file name error

    script:

    """
    # Need to decompress .gz file since fastq-pair doesn't support .gz files.

    mkdir match_pair_first_read
    cp $corrected_barcode_fastq match_pair_first_read/
    cp $read1_fastq match_pair_first_read/
    cd match_pair_first_read
    mv $corrected_barcode_fastq corrected_barcode.fq.gz
    mv $read1_fastq first_read_in_pair.fq.gz
    gzip -d corrected_barcode.fq.gz
    gzip -d first_read_in_pair.fq.gz
    fastq_pair $options.args corrected_barcode.fq first_read_in_pair.fq
    rm corrected_barcode.fq first_read_in_pair.fq
    gzip corrected_barcode.fq.paired.fq
    gzip first_read_in_pair.fq.paired.fq
    mv corrected_barcode.fq.paired.fq.gz first_read_matched_corrected_barcode.fq.gz

    cd ../

    mkdir match_pair_second_read
    cp $corrected_barcode_fastq match_pair_second_read/
    cp $read1_fastq match_pair_second_read/
    cd match_pair_second_read
    mv $corrected_barcode_fastq corrected_barcode.fq.gz
    mv $read1_fastq second_read_in_pair.fq.gz
    gzip -d corrected_barcode.fq.gz
    gzip -d second_read_in_pair.fq.gz
    fastq_pair $options.args corrected_barcode.fq second_read_in_pair.fq
    rm corrected_barcode.fq second_read_in_pair.fq
    gzip corrected_barcode.fq.paired.fq
    gzip second_read_in_pair.fq.paired.fq
    mv corrected_barcode.fq.paired.fq.gz second_read_matched_corrected_barcode.fq.gz

    """
}
