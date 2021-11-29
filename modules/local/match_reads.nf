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
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(corrected_barcode_fastq)

    output:
    tuple val(sample_name), path("match_pair_first_read/paired_R1_*.fastq.gz"), path("match_pair_second_read/paired_R2_*.fastq.gz"), path("match_pair_first_read/paired_barcode_*.fastq.gz"), path("match_pair_second_read/paired_barcode_*.fastq.gz"), emit: reads_2
    val sample_name, emit: sample_name
    // Be carefule of the duplicated staged file name error

    script:

    """
    # Need to decompress .gz file since fastq-pair doesn't support .gz files.

    mkdir match_pair_first_read
    cp $corrected_barcode_fastq match_pair_first_read/
    cp $read1_fastq match_pair_first_read/
    cd match_pair_first_read
    mv $corrected_barcode_fastq ${sample_name}_corrected_barcode.fq.gz
    mv $read1_fastq ${sample_name}_first_read_in_pair.fq.gz
    gzip -d ${sample_name}_corrected_barcode.fq.gz
    gzip -d ${sample_name}_first_read_in_pair.fq.gz
    fastq_pair $options.args ${sample_name}_corrected_barcode.fq ${sample_name}_first_read_in_pair.fq
    rm ${sample_name}_corrected_barcode.fq ${sample_name}_first_read_in_pair.fq
    gzip ${sample_name}_corrected_barcode.fq.paired.fq
    gzip ${sample_name}_first_read_in_pair.fq.paired.fq
    # mv ${sample_name}_corrected_barcode.fq.paired.fq.gz ${sample_name}_first_read_matched_corrected_barcode.fq.gz
    mv ${sample_name}_corrected_barcode.fq.paired.fq.gz paired_${barcode_fastq}
    mv ${sample_name}_first_read_in_pair.fq.paired.fq.gz paired_${read1_fastq}


    cd ../

    mkdir match_pair_second_read
    cp $corrected_barcode_fastq match_pair_second_read/
    cp $read2_fastq match_pair_second_read/
    cd match_pair_second_read
    mv $corrected_barcode_fastq ${sample_name}_corrected_barcode.fq.gz
    mv $read2_fastq ${sample_name}_second_read_in_pair.fq.gz
    gzip -d ${sample_name}_corrected_barcode.fq.gz
    gzip -d ${sample_name}_second_read_in_pair.fq.gz
    fastq_pair $options.args ${sample_name}_corrected_barcode.fq ${sample_name}_second_read_in_pair.fq
    rm ${sample_name}_corrected_barcode.fq ${sample_name}_second_read_in_pair.fq
    gzip ${sample_name}_corrected_barcode.fq.paired.fq
    gzip ${sample_name}_second_read_in_pair.fq.paired.fq
    # mv ${sample_name}_corrected_barcode.fq.paired.fq.gz ${sample_name}_second_read_matched_corrected_barcode.fq.gz
    mv ${sample_name}_corrected_barcode.fq.paired.fq.gz paired_${barcode_fastq}
    mv ${sample_name}_second_read_in_pair.fq.paired.fq.gz paired_${read2_fastq}

    """
}
