// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process DOWNLOAD_FROM_ENSEMBL {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_from_ensembl', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/miniconda3_xenial:0.1"

    // cache false

    input:
    val genome_name
    path dict_json

    output:
    path "*.fa.gz", emit: genome_fasta
    path "CHECKSUMS", emit: genome_md5
    val genome_name, emit: genome_name

    script:
    """
    md5_link=\$(get_download_url.py $dict_json $genome_name genome_md5sum)
    genome_link=\$(get_download_url.py $dict_json $genome_name genome)

    wget \$md5_link -o logfile.md5.txt
    wget \$genome_link -o logfile.genome.txt

    (cat \$(basename \$md5_link) | grep \$( basename \$genome_link) || true) > md5_to_check.txt

    if [ -s md5_to_check.txt ]
    then
      real=\$(cat md5_to_check.txt | cut -f 1,2 -d " ")
      measure=\$(sum \$( basename \$genome_link) | cut -f 1,2 -d " ")

      real1=\$(cat md5_to_check.txt | cut -f 1 -d " " | sed 's/^0*//')
      real2=\$(cat md5_to_check.txt | cut -f 2 -d " " | sed 's/^0*//')
      measure1=\$(sum \$( basename \$genome_link) | cut -f 1 -d " " | sed 's/^0*//')
      measure2=\$(sum \$( basename \$genome_link) | cut -f 2 -d " " | sed 's/^0*//')

      #if [ "\$real" == "\$measure" ]
      if [ "\$real1" == "\$measure1" ] && [ "\$real2" == "\$measure2" ]
      then
        exit 0
      else
        echo "\$real1,\$measure1,\$real2,\$measure2, md5check not pass"
        exit 1
      fi
    fi

    """
}
