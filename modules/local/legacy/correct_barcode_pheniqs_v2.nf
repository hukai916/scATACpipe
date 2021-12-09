// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CORRECT_BARCODE_PHENIQS {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode_pheniqs', publish_id:'') }
    container "hukai916/pheniqs_xenial:0.2"

    input:
    val sample_name
    path read1_fastq
    path read2_fastq
    path barcode_fastq
    path valid_barcode_frequency

    output:
    tuple val(sample_name), path("R1_*pheniqs*.fastq.gz"), path("R2_*pheniqs*.fastq.gz"), emit: reads_0
    path "summary_*.txt", emit: corrected_barcode_summary
    // barcode is not needed for pheniqs since corrected barcodes are added by default.

    script:

    """
    # step1, interleave read and index files
    pheniqs mux -R log_interleave.txt -i $read1_fastq -i $barcode_fastq -i $read2_fastq --output ${sample_name}.cram

    # step2, make a json config file
    barcode_length=\$((zcat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    make_json.py ${sample_name}_valid_barcode_frequency.txt ${sample_name}.cram 3 0::,2:: 1::\$barcode_length ${sample_name}.json

    # step3, run pheniqs
    pheniqs mux -R log_decode.txt --threads $task.cpus --decoding-threads $task.cpus --htslib-threads $task.cpus --config ${sample_name}.json --output ${sample_name}.corrected.bam

    # step4, extract fastq from pheniqs output bam
    samtools index ${sample_name}.corrected.bam
    bam2fastq.py ${sample_name}.corrected.bam corrected_${barcode_fastq}

    # rename output:
    sample_name=$read1_fastq
    outname="\${sample_name%%.*}"
    mv corrected_${barcode_fastq}*first_read_in_pair.fastq.gz \${outname}.pheniqs.fastq.gz
    sample_name=$read2_fastq
    outname="\${sample_name%%.*}"
    mv corrected_${barcode_fastq}*second_read_in_pair.fastq.gz \${outname}.pheniqs.fastq.gz

    # Note the noise param is determined by sequencer, can't be estimated; confidence: 0.99 (posterior possiblity), one run of pheniqs is okay to estimate the priors (since the invalid barcode are rare, this iteration is not a must.)


    """
}
