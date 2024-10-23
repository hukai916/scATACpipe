// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLIT_BED {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'split_bed', publish_id:'') }
    container "hukai916/miniconda3_v24.7_amd64_bio:0.1"

    input:
    tuple val(sample_name), path(fragment), path(tsv)

    output:
    path "split_*", emit: split_bed

    script:

    """
    tsv=($tsv)

    for (( i=0; i<\${#tsv[@]}; i++ )); do
      split_bed.py \${tsv[\$i]} $fragment
    done

    """
}
