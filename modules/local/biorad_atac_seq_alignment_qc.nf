// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process BIORAD_ATAC_SEQ_ALIGNMENT_QC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'biorad_alignment_qc', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/biorad_atac_seq_alignment_qc:0.1"
    // cache false

    input:
    val sample_name
    path alignments
    path bwa_fasta

    output:
    path "alignment_qc", emit: alignment_qc
    val sample_name, emit: sample_name

    script:

    """
    /runAlignmentQC.sh -i $alignments -r $bwa_fasta -o alignment_qc

    """
}
