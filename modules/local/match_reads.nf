// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process MATCH_READS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'match_reads', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/seqkit_0.16.1:0.1"

    // cache false

    input:
    val sample_name
    path corrected_barcode_fastq
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name
    // path "R2/barcode_corrected*fastq.gz", emit: barcode_fastq
    path corrected_barcode_fastq, emit: barcode_fastq
    path "R1/*.fastq.gz", emit: read1_fastq
    path "R2/*.fastq.gz", emit: read2_fastq

    script:

    """
    seqkit pair $options.args -1 $corrected_barcode_fastq -2 $read1_fastq -O R1
    rm R1/$corrected_barcode_fastq
    seqkit pair $options.args -1 R1/*.fastq.gz -2 $read2_fastq -O R2
    mv R2/$read1_fastq R1/
    # note that though the orders of R1 and R2 matches, the barcode may not.

    """
}
