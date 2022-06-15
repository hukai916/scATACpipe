// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TAG_BAM {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'tag_bam', publish_id:'') }

    container "hukai916/sinto_xenial:0.2"

    input:
    tuple val(sample_name), val(chunk_name), path(tagfile), path(bam)
    // val sample_name
    // path tagfile
    // path bam

    output:
    val sample_name, emit: sample_name
    path "*.tag.bam", emit: bam

    script:

    """
    # v1: use more memory
    # samtools view -h $bam | awk 'BEGIN{FS=OFS="\t"} NR == FNR { if (\$3 != "undetermined") {tag_dict[\$1] = \$3}; next} { if (/^@/) {print} else if (\$1 in tag_dict) {print \$0, "CB:Z:" tag_dict[\$1]}}' $tagfile - | samtools view -b -h - -o ${chunk_name}.tag.bam

    # v2: use less memory but assuming that all qnames in the BAM files are in the tagfile, which is not true since tagfile only contains valid barcodes.
    # samtools view -h $bam | awk 'BEGIN{FS=OFS="\t"} NR == FNR { if (\$2 != \$3) { if (\$3 == "undetermined") { tag_dict[\$1] = 0 } else { tag_dict[\$1] = \$3 } }; next} { if (/^@/) {print} else { split(\$1, q_name, ":"); print \$0, "CB:Z:" q_name[1]}}'  $tagfile - | samtools view -h -b - -o ${chunk_name}.tag.bam

    # v3: use less memory and without assumption in v2
    samtools view -h $bam | awk 'BEGIN{FS=OFS="\t"} NR == FNR { if (\$3 != "undetermined") { if (\$2 == \$3) { tag_dict[\$1] = 1 } else { tag_dict[\$1] = \$3 } }; next} { if (/^@/) {print} else if (\$1 in tag_dict) { if (tag_dict[\$1] != 1) { print \$0, "CB:Z:" tag_dict[\$1] } else { split(\$1, q_name, ":"); print \$0, "CB:Z:" q_name[1] } } }'  $tagfile - | samtools view -h -b - -o ${chunk_name}.tag.bam

    """
}
