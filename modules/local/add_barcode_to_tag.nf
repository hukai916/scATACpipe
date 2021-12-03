// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADD_BARCODE_TO_TAG {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'add_barcode_to_tag', publish_id:'') }
    container "hukai916/sinto_xenial:0.2"

    input:
    val sample_name
    path bam

    output:
    val sample_name, emit: sample_name
    path "*.barcode_tagged.bam", emit: bam

    script:

    """
    # Copy cell barcode from readname to tag CB (assuming barcode matches /[^:]*/):
    samtools view -h $bam | awk 'BEGIN { OFS = "\t"} match(\$1, /[^:]*/) { print \$0, "CR:Z:"substr(\$1, RSTART, RLENGTH) }' | samtools view -o ${bam.baseName}.barcode_tagged.bam

    """

}
