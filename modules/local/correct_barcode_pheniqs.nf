// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CORRECT_BARCODE_PHENIQS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode_pheniqs', publish_id:'') }
    container "hukai916/pheniqs_xenial:0.1"

    input:
    val sample_name
    path barcode_fastq
    path barcode_whitelist
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name
    path "summary_*.txt", emit: corrected_barcode_summary
    path "barcode_corrected*.R1.fastq.gz", emit: read1_fastq
    path "barcode_corrected*.R2.fastq.gz", emit: read2_fastq
    path "barcode_corrected*.R3.fastq.gz", emit: corrected_barcode

    script:

    """

    # step1, interleave read and index files
    pheniqs mux -R log_interleave.txt -i $read1_fastq -i $barcode_fastq -i $read2_fastq --output ${sample_name}.cram

    # step2, retrieve valid barcode pool and concentration in raw counts
    get_barcode_pool.py $barcode_whitelist $barcode_fastq $options.read_count_cutoff valid_barcode_pool.txt

    # step3, make a json config file
    barcode_length=\$(awk 'NR==1 {print length(\$(NF-1)); exit}' valid_barcode_pool.txt)
    make_json.py valid_barcode_pool.txt ${sample_name}.cram 3 0::,2:: 1::\$barcode_length ${sample_name}.json

    # step4, run pheniqs
    pheniqs mux -R log_decode.txt --threads $task.cpus --decoding-threads $task.cpus --htslib-threads $task.cpus --config ${sample_name}.json --output ${sample_name}.bam

    # step5, extract fastq from pheniqs output bam
    bam2fastq.py ${sample_name}.bam barcode_corrected_${sample_name}

    # Note the noise param is determined by sequencer, can't be estimated; confidence: 0.99 (posterior possiblity), one run of pheniqs is okay to estimate the priors (since the invalid barcode are rare, this iteration is not a must.)


    """
}
