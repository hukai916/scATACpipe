// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWA_INDEX {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'bwa_index', publish_id:'') }
    container "hukai916/bwa_xenial:0.1"

    input:
    path genome_fasta

    output:
    path "bwa_index", emit: bwa_index_folder

    script:
    genome_basename = genome_fasta.getName()

    """
    mkdir bwa_index
    ln $genome_fasta bwa_index/
    bwa index $options.args -a bwtsw bwa_index/$genome_basename

    """
}
