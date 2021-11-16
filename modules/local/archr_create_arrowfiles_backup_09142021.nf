// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_CREATE_ARROWFILES {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_create_arrowfiles', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc:0.5"

    // cache false

    input:
    tuple val(sample_name), path(fragment)
    // val sample_name
    // path fragment
    val archr_genome
    val archr_thread

    output:
    val sample_name, emit: sample_name
    path "QualityControl_*", emit: quality_control
    path "*.arrow", emit: arrowfile

    script:
    // for unknown reason, #!/usr/bin/R + direct R codes won't work
    """
    echo '
    library(ArchR)
    addArchRGenome("$archr_genome")
    addArchRThreads(threads = $archr_thread)

    inputFiles <- "$fragment"
    names(inputFiles) <- "$sample_name"

    ArrowFiles <- createArrowFiles(
      inputFiles = inputFiles,
      sampleNames = names(inputFiles),
      threads = 1,
      QCDir = paste0("QualityControl_", "$sample_name"),
      subThreading = FALSE,
      $options.args
    )' > run.R

    Rscript run.R

    # Convert to jpeg:
    mkdir QualityControl_$sample_name/jpeg
    x=( \$(find ./QualityControl_$sample_name -name "*.pdf") )
    for item in "\${x[@]}"
    do
      filename=\$(basename -- "\$item")
      filename="\${filename%.*}"
      pdftoppm -jpeg -r 300 \$item ./QualityControl_$sample_name/jpeg/\$filename
      convert -append ./QualityControl_$sample_name/jpeg/\${filename}* ./QualityControl_$sample_name/jpeg/\${filename}.jpg
      rm ./QualityControl_$sample_name/jpeg/\${filename}-*.jpg
    done

    """
}
