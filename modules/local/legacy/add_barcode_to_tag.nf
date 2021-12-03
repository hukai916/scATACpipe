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
    # Copy cell barcode from readname to tag CB:
    # Prepare tag file:
    samtools view $bam | awk 'BEGIN { OFS = "\t"} match(\$1, /[^:]*/) { print substr(\$1, RSTART, RLENGTH), "CB", substr(\$1, RSTART, RLENGTH)}' > ${bam.baseName}_tag.tsv

    # Add barcode to CB tag:
    samtools index $bam
    sinto addtags -p $task.cpus -m tag -b $bam -f ${bam.baseName}_tag.tsv -o ${bam.baseName}.barcode_tagged.bam

    """

}
