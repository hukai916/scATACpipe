// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_PEAK2GENELINKAGE_CLUSTERS2 {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_peak2genelinkage_clusters2', publish_id:'') }

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
    path "Plots/Heatmap-Marker-Genes-with-Peak2GeneLinks.pdf", emit: heatmap_marker_genes_with_peaks2genelinks
    path "Plots/Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", emit: plot_tracks_marker_genes_with_peak2genelinks
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "archr_peak2genelinkage_clusters2", emit: res_dir
    path "report_jpeg/archr_peak2genelinkage_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- addPeak2GeneLinks(
      ArchRProj = proj,
      reducedDims = "IterativeLSI"
      )

    p2g <- getPeak2GeneLinks(
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
      loops = getPeak2GeneLinks(proj2)
      )
    plotPDF(plotList = p,
      name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",
      ArchRProj = NULL,
      addDOC = FALSE, width = 5, height = 5)

    p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Clusters2")
    plotPDF(p,
      name = "Heatmap-Marker-Genes-with-Peak2GeneLinks.pdf",
      ArchRProj = NULL,
      addDOC = FALSE, width = 5, height = 12)

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

    # Copy to res_dir:
    mkdir archr_peak2genelinkage_clusters2
    cp -r Plots archr_peak2genelinkage_clusters2/

    # For reporting:
    mkdir -p report_jpeg/archr_peak2genelinkage_clusters2
    cp -r archr_peak2genelinkage_clusters2 report_jpeg/archr_peak2genelinkage_clusters2

    """
}
