// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_MARKER_GENE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_marker_gene', publish_id:'') }

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
    path "proj_marker_gene.rds", emit: archr_project
    path "Plots/GeneScores-Marker-Heatmap.pdf", emit: pdf_genescores_marker_heatmap
    path "Plots/Plot-UMAP-Marker-Genes-WO-Imputation.pdf", emit: pdf_umap_markder_genes_wo_imputation
    path "Plots/Plot-UMAP-Marker-Genes-W-Imputation.pdf", emit: pdf_umap_markder_genes_w_imputation
    path "Plots/Plot-Tracks-Marker-Genes.pdf", emit: pdf_tracks_marker_genes
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_marker_gene", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    # Find marker genes:
    markersGS <- getMarkerFeatures(
      ArchRProj = proj,
      $options.args
    )

    markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

    # Draw heatmap:
    markerGenes  <- c($options.marker_genes)
    heatmapGS <- markerHeatmap(
      seMarker = markersGS,
      cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
      labelMarkers = markerGenes,
      transpose = TRUE
    )
    plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    # Plot marker genes on embedding:
    p <- plotEmbedding(
      ArchRProj = proj,
      name = markerGenes,
      imputeWeights = NULL,
      $options.args2
    )
    plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    # Add Marker Genes Imputation and replot:
    proj2 <- addImputeWeights(proj)
    saveRDS(proj2, file = "proj_marker_gene.rds")

    p <- plotEmbedding(
      ArchRProj = proj2,
      name = markerGenes,
      imputeWeights = getImputeWeights(proj2),
      $options.args2
    )
    plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    # Plot: track plotting with ArchRBrowser
    p <- plotBrowserTrack(
      ArchRProj = proj2,
      geneSymbol = markerGenes,
      $options.args3
    )
    plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

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
    mkdir -p report_jpeg/archr_marker_gene
    cp -r Plots/jpeg report_jpeg/archr_marker_gene

    """
}
