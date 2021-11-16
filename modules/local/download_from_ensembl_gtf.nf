// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process DOWNLOAD_FROM_ENSEMBL_GTF {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'download_from_ensembl_gtf', publish_id:'') }

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
    path "*.gtf.gz", emit: gtf
    path "CHECKSUMS", emit: genome_md5
    val genome_name, emit: genome_name

    script:
    """
    md5_link=\$(get_download_url.py $dict_json $genome_name gtf_md5sum)
    gtf_link=\$(get_download_url.py $dict_json $genome_name gtf)

    wget \$md5_link -o logfile.md5.txt
    wget \$gtf_link -o logfile.gtf.txt

    (cat \$(basename \$md5_link) | grep \$( basename \$gtf_link) || true) > md5_to_check.txt

    if [ -s md5_to_check.txt ]
    then # note that also need to remove leading 0s if any
      real1=\$(cat md5_to_check.txt | cut -f 1 -d " " | sed 's/^0*//')
      real2=\$(cat md5_to_check.txt | cut -f 2 -d " " | sed 's/^0*//')
      measure1=\$(sum \$( basename \$gtf_link) | cut -f 1 -d " " | sed 's/^0*//')
      measure2=\$(sum \$( basename \$gtf_link) | cut -f 2 -d " " | sed 's/^0*//')

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
