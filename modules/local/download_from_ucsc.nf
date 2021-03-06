// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DOWNLOAD_FROM_UCSC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_from_ucsc', publish_id:'') }
    container "hukai916/miniconda3_xenial:0.1"

    input:
    val genome_name
    path dict_json

    output:
    path "*.fa.gz", emit: genome_fasta
    path "md5sum.txt", emit: genome_md5
    val genome_name, emit: genome_name

    script:

    """
    md5_link=\$(get_download_url.py $dict_json $genome_name md5sum)
    genome_link=\$(get_download_url.py $dict_json $genome_name genome)

    # if [[ md5_link -eq 0 ]]; then { echo "Genome not supported!"; exit; } fi

    wget --no-check-certificate \$md5_link -o logfile.md5.txt
    wget --no-check-certificate \$genome_link -o logfile.genome.txt

    (cat \$(basename \$md5_link) | grep \$( basename \$genome_link) || true) > md5_to_check.txt

    if [ -s md5_to_check.txt ]
    then
      md5sum -c md5_to_check.txt
    fi

    """
}
