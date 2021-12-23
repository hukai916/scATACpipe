// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process QUALIMAP {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'qualimap', publish_id:'') }
    container "hukai916/qualimap_xenial:0.1"

    input:
    tuple val(sample_name), path(bam)

    output:
    val sample_name, emit: sample_name
    path "$sample_name", emit: bamqc

    script:

    """
    mem=\$(echo '$task.memory' | grep -o -E '[0-9]+')
    qualimap bamqc $options.args --java-mem-size=\${mem}G -bam $bam -outdir ${sample_name}

    """
}
