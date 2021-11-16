// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process CORRECT_BARCODE_PHENIQS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode_pheniqs', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/pheniqs_xenial:0.1"

    // cache false

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
    barcode_length=\$(awk 'FNR==1 {print length(\$(NF-1))}' valid_barcode_pool.txt)
    make_json.py valid_barcode_pool.txt ${sample_name}.cram 3 0::,2:: 1::\$barcode_length ${sample_name}.json

    # step4, run pheniqs
    pheniqs mux -R log_decode.txt --threads $task.cpus --decoding-threads $task.cpus --htslib-threads $task.cpus --config ${sample_name}.json --output ${sample_name}.bam

    # step5, extract fastq from pheniqs output bam
    bam2fastq.py ${sample_name}.bam barcode_corrected_${sample_name}

    # step1-1 calculate barcode pool and frequency using our own script.

    # step2 uses step1-1. noise: 0.01 (determined by sequencer, can't be estimated), confidence: 0.99 (posterior possiblity), one run of pheniqs is okay to estimate the priors (since the invalid barcode are rare, this iteration is not a must.)

    # step2, prepare a json config file for pheniqs
    # use valid bacode to calculate frequency and pool.
    # share with Haibo a demo command to run 10xgenomics.

    """
}
