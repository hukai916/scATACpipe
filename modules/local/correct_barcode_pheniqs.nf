// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CORRECT_BARCODE_PHENIQS {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'correct_barcode_pheniqs', publish_id:'') }
    container "hukai916/pheniqs_xenial:0.2"

    input:
    tuple val(sample_name), path(barcode_fastq), path(valid_barcode_counts)

    output:
    tuple val(sample_name), val(chunk_name), path("*.tag.tsv"), emit: sample_name_chunk_name_tagfile
    path "summary_*.txt", emit: corrected_barcode_summary

    script:
    // In split_fastq is not called, chunk_name remains the same as sample_name.
    chunk_name = barcode_fastq.name.split("\\.")[0..-3].join(".").split("barcode_").join() // get rid of suffix ".fastq.gz", then remove leading "barcode_" if there.

    """
    # Step1, interleave read and index files:
    cp $barcode_fastq temp.fastq.gz # pheniqs do not accept duplicate filenames
    pheniqs mux -R log_interleave.txt -i temp.fastq.gz -i $barcode_fastq --output ${chunk_name}.cram

    # Step2, make a json config file (use a minial 0:0:1 as output R1 to save I/O):
    barcode_length=\$((zcat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    make_json.py $valid_barcode_counts ${chunk_name}.cram 2 0:0:1 1::\$barcode_length ${chunk_name}.json

    # Step3, run pheniqs:
    pheniqs mux -R log_decode.txt --threads $task.cpus --decoding-threads $task.cpus --htslib-threads $task.cpus --config ${chunk_name}.json --output ${chunk_name}.corrected.bam

    # Step4, extract a tag file:
    samtools index ${chunk_name}.corrected.bam
    extract_tag.py ${chunk_name}.corrected.bam BC,RG ${chunk_name}.tag.tsv
    # cat ${chunk_name}.tag.tsv | awk 'BEGIN { OFS = "\\t"} { print \$1,"CB",\$2 }' > ${chunk_name}.tagfile_sinto.tsv

    # Step5, print stats:
    valid_read_num=\$(cat ${chunk_name}.tag.tsv | awk '{ if (\$2 == \$3) print \$0 }' | wc -l)
    discard_read_num=\$(cat ${chunk_name}.tag.tsv | awk '{ if (\$3 == "undetermined") print \$0 }' | wc -l)
    rescued_read_num=\$(cat ${chunk_name}.tag.tsv | grep -v "undetermined" | awk '{ if (\$2 != \$3) print \$0 }' | wc -l)

    echo "Summary (correct_barcode): total valid: "\${valid_read_num}"; total corrected: "\${rescued_read_num}"; total discarded: "\${discard_read_num}"." > summary_${chunk_name}.txt

    """
}
