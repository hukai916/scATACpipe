// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREP_BAM {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'prep_bam', publish_id:'') }
    container "hukai916/sinto_kai:0.1"

    input:
    val sample_name
    path bam

    output:
    val sample_name, emit: sample_name
    path "*.prep.bam", emit: bam
    path "summary_prep_bam_*.txt", emit: summary_prep_bam

    script:

    """
    # Prepare bam files for combine_bam.nf module:
    # 1. Copy cell barcode from name to tag (default to CR (cell barcode reported by sequencer));
    # 2. Shift reads to compensate for Tn5;
    # 3. Extend soft clips;
    # 4. Add RG tag using barcode.

    # For 1-3:
    prep_bam.py --barcode_tag_to_add CR --inbam $bam --outbam ${sample_name}.prep.temp.bam $options.args

    # For 4:
    samtools view ${sample_name}.prep.temp.bam | awk 'BEGIN { OFS = "\t"} match(\$0, /CR:Z:[a-zA-Z]*/) { print substr(\$0, RSTART+5, RLENGTH-5)}' | sort | uniq > ${sample_name}.uniq.CR.txt
    samtools index ${sample_name}.prep.temp.bam
    sinto tagtorg -b ${sample_name}.prep.temp.bam --tag CR -f ${sample_name}.uniq.CR.txt -o ${sample_name}.prep.bam

    # Clean up:
    rm ${sample_name}.prep.temp.bam ${sample_name}.prep.temp.bam.bai

    """

}
