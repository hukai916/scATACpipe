// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHROMAP_INDEX {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'chromap_index', publish_id:'') }
    container "hukai916/chromap_xenial:0.1"

    input:
    tuple path(genome_fasta), val(genome_name)

    output:
    path "chromap_index_*", emit: index_file

    script:

    """
    chromap $options.args -i -r $genome_fasta -o chromap_index_$genome_name

    """
}
