// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_VALID_BARCODE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_valid_barcode', publish_id:'') }
    container "hukai916/umitools_xenial:0.1"

    input:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)
    path whitelist_barcode

    output:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq), emit: reads
    path "*valid_barcode_frequency.txt", emit: valid_barcode_frequency
    path "*valid_barcode.txt", emit: valid_barcode

    script:

    """
    barcode_length=\$((zcat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    printf -v bc_pattern '%0.sC' \$(seq 1 \$barcode_length)
    umi_tools whitelist $options.args -I $barcode_fastq --bc-pattern \$bc_pattern --error-correct-threshold 0 | grep -v "#" > valid_barcode_frequency_raw.txt

    if [[ $whitelist_barcode == file_token.txt ]]; then
      cat valid_barcode_frequency_raw.txt | cut -f 1 > ${sample_name}_valid_barcode.txt
    else
      get_valid_barcode.py valid_barcode_frequency_raw.txt $whitelist_barcode ${sample_name}_valid_barcode.txt ${sample_name}_valid_barcode_frequency.txt
    fi

    """
}
