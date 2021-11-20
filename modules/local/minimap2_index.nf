// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIMAP2_INDEX {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'minimap2_index', publish_id:'') }
    container "hukai916/minimap2_xenial:0.1"

    input:
    path genome_fasta

    output:
    path "*.mmi", emit: minimap2_index

    script:
    genome_basename = genome_fasta.getName()

    """
    minimap2 -t $task.cpus $options.args -d ${genome_basename}.mmi $genome_fasta

    """
}
