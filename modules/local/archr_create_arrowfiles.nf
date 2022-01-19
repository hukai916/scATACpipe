// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_CREATE_ARROWFILES {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_create_arrowfiles', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    tuple val(sample_name), path(fragment)
    val archr_genome
    val archr_thread

    output:
    val sample_name, emit: sample_name
    path "QualityControl_*", emit: quality_control
    path "*.arrow", emit: arrowfile
    path "report_*", emit: report

    beforeScript params.archr_beforescript

    script:

    """
    echo \$HDF5_USE_FILE_LOCKING > test.txt

    echo '
    library(ArchR)

    inputFiles <- "$fragment"
    names(inputFiles) <- "$sample_name"

    addArchRGenome("$archr_genome")
    addArchRThreads(threads = $archr_thread)

    ArrowFiles <- createArrowFiles(
      inputFiles = inputFiles,
      sampleNames = names(inputFiles),
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
    mkdir -p report_archr_create_arrowfiles_$sample_name/archr_create_arrowfiles
    cp -r QualityControl_$sample_name/jpeg report_archr_create_arrowfiles_$sample_name/archr_create_arrowfiles

    """
}
