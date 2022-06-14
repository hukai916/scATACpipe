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

    wget --user-agent=" Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36" \$md5_link -o logfile.md5.txt
    wget --user-agent=" Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36" \$genome_link -o logfile.genome.txt

    (cat \$(basename \$md5_link) | grep \$( basename \$genome_link) || true) > md5_to_check.txt

    if [ -s md5_to_check.txt ]
    then
      md5sum -c md5_to_check.txt
    fi

    """
}
