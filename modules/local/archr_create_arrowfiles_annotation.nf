// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_CREATE_ARROWFILES_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_create_arrowfiles_annotation', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_archr:0.1"

    // cache false

    input:
    tuple val(sample_name), path(fragment)
    path gene_annotation
    path genome_annotation
    path user_rlib
    val archr_thread

    output:
    val sample_name, emit: sample_name
    path "QualityControl_*", emit: quality_control
    path "*.arrow", emit: arrowfile
    path "report_*", emit: report

    script:
    // for unknown reason, #!/usr/bin/R + direct R codes won't work

    """
    echo '
    library(ArchR)

    # Include the installed custom BSgenome if supplied:
    if (!("$user_rlib" == "file_token.txt")) {
      .libPaths("user_rlib")
    }

    inputFiles <- "$fragment"
    names(inputFiles) <- "$sample_name"

    genomeAnnotation <- readRDS("$genome_annotation")
    geneAnnotation <- readRDS("$gene_annotation")

    ArrowFiles <- createArrowFiles(
      inputFiles = inputFiles,
      sampleNames = names(inputFiles),
      geneAnnotation = geneAnnotation,
      genomeAnnotation = genomeAnnotation,
      threads = 1,
      QCDir = paste0("QualityControl_", "$sample_name"),
      subThreading = FALSE,
      $options.args
    )
' > run.R

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

    # For reporting:
    mkdir -p report_archr_create_arrowfiles_annotation_$sample_name/archr_create_arrowfiles_annotation
    cp -r QualityControl_$sample_name/jpeg report_archr_create_arrowfiles_annotation_$sample_name/archr_create_arrowfiles_annotation


    """
}
