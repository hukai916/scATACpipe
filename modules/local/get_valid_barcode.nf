// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_VALID_BARCODE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_valid_barcode', publish_id:'') }
    container "hukai916/r_util:0.2"

    input:
    val sample_name
    path barcode_fastq
    path whitelist_barcode

    output:
    val sample_name, emit: sample_name
    path barcode_fastq, emit: barcode_fastq
    path "*_valid_barcode_counts_fastq.txt", emit: valid_barcode_counts_fastq
    path "*_valid_barcodes_dedup_bam.txt", emit: valid_barcodes

    script:

    """
    # Generate three outfiles:
    #  Outfile1: Barcode read counts from dedup bam: -> inflection point -> valid barcode counts (outfile2).
    #  Outfile3: Barcode counts from raw barcode fastq: -> raw barcode counts -> pheniqs barcode correction.

    # For outfile1:
    samtools view ${sample_name}.dedup.bam | awk 'BEGIN { OFS = "\t" } match(\$1, /[^:]*/) { print substr(\$1, RSTART, RLENGTH) }' | sort | uniq -c | awk '{print \$2, \$1}' > ${sample_name}_barcode_counts_dedup_bam.txt

    # For outfile2:
    get_valid_barcode_inflection.R --freq ${sample_name}_barcode_counts_dedup_bam.txt --outfile ${sample_name}_valid_barcode_counts_dedup_bam_temp.txt

    if [[ $whitelist_barcode == file_token.txt ]]; then
      cat ${sample_name}_valid_barcode_counts_dedup_bam_temp.txt | cut -f 1 > ${sample_name}_valid_barcodes_dedup_bam.txt
      mv ${sample_name}_valid_barcode_counts_dedup_bam_temp.txt ${sample_name}_valid_barcode_counts_dedup_bam.txt
    else
      join -1 1 -2 1 <(sort $whitelist_barcode) <(sort ${sample_name}_valid_barcode_counts_dedup_bam_temp.txt) > ${sample_name}_valid_barcode_counts_dedup_bam.txt
      rm ${sample_name}_valid_barcode_counts_dedup_bam_temp.txt
      cat ${sample_name}_valid_barcode_counts_dedup_bam.txt | cut -f 1 > ${sample_name}_valid_barcodes_dedup_bam.txt
    fi

    # For outfile3:
    zcat $barcode_fastq | awk 'NR%4==2 {print}' | sort | uniq -c | awk '{print \$2, \$1}' > ${sample_name}_barcode_counts_fastq.txt
    join -1 1 -2 1 <(sort ${sample_name}_valid_barcodes_dedup_bam.txt) <(sort ${sample_name}_barcode_counts_fastq.txt) > ${sample_name}_valid_barcode_counts_fastq.txt

    """
}
