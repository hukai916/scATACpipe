// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_VALID_BARCODE_CELLRANGER {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_valid_barcode_cellranger', publish_id:'') }
    container "hukai916/umitools_xenial:0.1"

    input:
    tuple val(sample_name), path(sample_files), path(bam), path(fragment)
    path whitelist_barcode

    output:
    tuple val(sample_name), path(sample_files), path(bam), path(fragment), emit: sample
    path "valid_barcode.txt", emit: valid_barcode

    script:

    """
    cat *_R2_001.fastq.gz > barcode.fastq.gz

    barcode_length=\$((zcat barcode.fastq.gz || true) | awk 'NR==2 {print length(\$0); exit}')
    printf -v bc_pattern '%0.sC' \$(seq 1 \$barcode_length)
    umi_tools whitelist $options.args -I barcode.fastq.gz --bc-pattern \$bc_pattern --error-correct-threshold 0 | grep -v "#" > valid_barcode_frequency.txt

    if [[ $whitelist_barcode == file_token.txt ]]; then
      cat valid_barcode_frequency.txt | cut -f 1 > valid_barcode.txt
    else
      get_valid_barcode.py valid_barcode_frequency.txt $whitelist_barcode valid_barcode.txt
    fi

    """
}
