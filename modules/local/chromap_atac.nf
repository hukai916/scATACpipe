// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHROMAP_ATAC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'chromap_atac', publish_id:'') }
    container "hukai916/chromap_xenial:0.1"

    input:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq), path(whitelist_barcode)
    path genome_fasta
    path genome_index
    val use_whitelist

    output:
    tuple val(sample_name), path("chromap_fragment_*"), emit: sample_name_fragments
    path "chromap_fragment_*", emit: fragments

    script:

    """
    echo "test"


    """
}
