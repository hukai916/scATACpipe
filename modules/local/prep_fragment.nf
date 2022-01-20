// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREP_FRAGMENT {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'prep_fragment', publish_id:'') }
    container "hukai916/miniconda3_bio:0.1"

    input:
    tuple val(sample_name), path(fragment)
    path gtf

    output:
    tuple val(sample_name), path("final.*.fragment.tsv.gz"), emit: fragments

    script:

    """
    if [[ $fragment == *.gz ]]; then
      gunzip -c $fragment > fragment.bed
    else
      if [[ $fragment != fragment.bed ]]; then
        mv $fragment fragment.bed
      fi
    fi

    if [[ $gtf == *.gz ]]; then
      gunzip -c $gtf > annotation.gtf
    else
      if [[ $gtf != annotation.gtf ]]; then
        mv $gtf annotation.gtf
      fi
    fi

    # add 'chr' prefix to fragment bed to match PREP_GTF
    cat fragment.bed | awk 'BEGIN{FS=OFS="\\t"}/^#/{print; next} NF>=3{\$1 = gensub(/(chr)?(.+).*/, "chr\\\\2", "g", \$1); print}' > chrPrefixed.fragment.bed

    # make sure bed col1 is a subset of gtf col1:
    extract_bed.py chrPrefixed.fragment.bed annotation.gtf | bgzip > final.${sample_name}.fragment.tsv.gz

    """
}
