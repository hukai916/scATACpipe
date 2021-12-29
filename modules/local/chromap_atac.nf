// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHROMAP_ATAC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'chromap_atac', publish_id:'') }
    container "hukai916/chromap_xenial:0.1"

    input:
    tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq), path(whitelist_barcode)
    path genome_fasta
    path genome_index
    val use_whitelist

    // val sample_name
    // path sample_fastqs
    // path genome_fasta
    // path index
    // path whitelist_barcode

    output:
    path "chromap_fragment_*", emit: fragment

    script:

    """
    if [[ $use_whitelist == false ]]; then
      option_whitelist=''
    else
      if [[ $whitelist_barcode == *.gz ]]; then
        gunzip -c $whitelist_barcode > whitelist_${sample_name}.txt
      else
        mv $whitelist_barcode whitelist_${sample_name}.txt
      fi
      option_whitelist='--barcode-whitelist whitelist_${sample_name}.txt'
    fi

    chromap --preset atac \
    $options.args \
    -t $task.cpus \
    -x $genome_index \
    -r $genome_fasta \
    -1 $read1_fastq \
    -2 $read2_fastq \
    -o chromap_fragment_${sample_name}.bed \
    -b $barcode_fastq \
    \$option_whitelist

    """
}