// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_CREATE_ARROWFILES_ANNOTATION {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_create_arrowfiles_annotation', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

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
    path "report_jpeg/archr_create_arrowfiles_annotation_*", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

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
      QCDir = paste0("QualityControl_", "$sample_name"),
      subThreading = FALSE,
      $options.args
    )
    ' > run.R

    Rscript run.R

    # Convert to jpeg:
    mkdir -p QualityControl_$sample_name/jpeg
    x=( \$(find ./QualityControl_$sample_name -name "*.pdf") )
    for item in \${x[@]+"\${x[@]}"}
    do
      {
        filename=\$(basename -- "\$item")
        filename="\${filename%.*}"
        pdftoppm -jpeg -r 300 \$item ./QualityControl_$sample_name/jpeg/\$filename
        convert -append ./QualityControl_$sample_name/jpeg/\${filename}* ./QualityControl_$sample_name/jpeg/\${filename}.jpg
        rm ./QualityControl_$sample_name/jpeg/\${filename}-*.jpg
      } || {
        echo "Pdf to jpeg failed!" > bash.log
      }
    done

    # For reporting:
    mkdir -p ./report_jpeg/archr_create_arrowfiles_annotation_$sample_name
    cp -r ./QualityControl_$sample_name/jpeg report_jpeg/archr_create_arrowfiles_annotation_$sample_name || :

    mkdir ./report_jpeg/archr_create_arrowfiles_annotation_$sample_name/pdf
    cp ./QualityControl_$sample_name/*/*.pdf report_jpeg/archr_create_arrowfiles_annotation_$sample_name/pdf || :

    """
}
