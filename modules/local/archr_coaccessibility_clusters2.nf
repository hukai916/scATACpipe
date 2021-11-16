// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_COACCESSIBILITY_CLUSTERS2 {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_coaccessibility_clusters2', publish_id:'') }

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
    path archr_project

    output:
    path "archr_coaccessibility_project.rds", emit: archr_project
    path "Plots/Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", emit: plot_tracks_marker_genes_with_coaccessibility
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_coaccessibility_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- addCoAccessibility(
      ArchRProj = proj,
      reducedDims = "IterativeLSI"
      )
    saveRDS(proj2, file = "archr_coaccessibility_project.rds")

    cA <- getCoAccessibility(
      ArchRProj = proj2,
      returnLoops = TRUE,
      $options.args
      )

    markerGenes <- c($options.marker_genes)

    p <- plotBrowserTrack(
      ArchRProj = proj2,
      groupBy = "Clusters2",
      geneSymbol = markerGenes,
      upstream = 50000,
      downstream = 50000,
      loops = getCoAccessibility(proj2)
      )

      plotPDF(plotList = p,
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
        ArchRProj = NULL,
        addDOC = FALSE, width = 5, height = 5
        )

    ' > run.R

    Rscript run.R

    # Convert to jpeg:
    mkdir Plots/jpeg
    x=( \$(find ./Plots -name "*.pdf") )
    for item in "\${x[@]}"
    do
      filename=\$(basename -- "\$item")
      filename="\${filename%.*}"
      pdftoppm -jpeg -r 300 \$item ./Plots/jpeg/\$filename
      convert -append ./Plots/jpeg/\${filename}* ./Plots/jpeg/\${filename}.jpg
      rm ./Plots/jpeg/\${filename}-*.jpg
    done

    # For reporting:
    mkdir -p report_jpeg/archr_coaccessibility_clusters2
    cp -r Plots/jpeg report_jpeg/archr_coaccessibility_clusters2

    """
}
