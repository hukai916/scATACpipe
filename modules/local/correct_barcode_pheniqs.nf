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
    val sample_name
    // path read1_fastq
    // path read2_fastq
    path barcode_fastq
    path valid_barcode_counts

    output:
    // tuple val(sample_name), path("R1_*pheniqs*.fastq.gz"), path("R2_*pheniqs*.fastq.gz"), emit: reads_0
    val sample_name, emit: sample_name
    path "*.tag.tsv", emit: tagfile
    path "summary_*.txt", emit: corrected_barcode_summary

    script:

    """
    # Step1, interleave read and index files:
    cp $barcode_fastq temp.fastq.gz # pheniqs do not accept duplicate filenames
    pheniqs mux -R log_interleave.txt -i temp.fastq.gz -i $barcode_fastq --output ${sample_name}.cram

    # Step2, make a json config file (use a minial 0:0:1 as output R1 to save I/O):
    barcode_length=\$((zcat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    make_json.py $valid_barcode_counts ${sample_name}.cram 2 0:0:1 1::\$barcode_length ${sample_name}.json

    # Step3, run pheniqs:
    pheniqs mux -R log_decode.txt --threads $task.cpus --decoding-threads $task.cpus --htslib-threads $task.cpus --config ${sample_name}.json --output ${sample_name}.corrected.bam

    # Step4, extract a tag file:
    samtools index ${sample_name}.corrected.bam
    extract_tag.py ${sample_name}.corrected.bam BC,RG ${sample_name}.tag.tsv
    # cat ${sample_name}.tag.txt | awk 'BEGIN { OFS = "\\t"} { print \$1,"CB",\$2 }' > ${sample_name}.tagfile_sinto.tsv

    # Step5, print stats:
    valid_read_num=\$(cat ${sample_name}.tag.txt | awk '{ if (\$1 == \$2) print \$0 }' | wc -l)
    discard_read_num=\$(cat ${sample_name}.tag.txt | awk '{ if (\$2 == "undetermined") print \$0 }' | wc -l)
    rescued_read_num=\$(cat ${sample_name}.tag.txt | grep -v "undetermined" | awk '{ if (\$1 != \$2) print \$0 }' | wc -l)

    echo "Summary (correct_barcode): total valid: "\${valid_read_num}"; total corrected: "\${rescued_read_num}"; total discarded: "\${discard_read_num}"." > summary_${sample_name}.txt

    """
}
