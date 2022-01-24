// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_MARKER_GENE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_marker_gene', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val archr_thread

    output:
    path "proj_marker_gene.rds", emit: archr_project
    path "marker_list.txt", emit: marker_list
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_marker_gene", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(stringr)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    # Find marker genes: default to use Seurat
    if ("Harmony" %in% names(proj@reducedDims)) {
      markersGS <- getMarkerFeatures(
        ArchRProj = proj,
        groupBy = "Clusters_Seurat_Harmony",
        $options.args
      )
    } else {
      markersGS <- getMarkerFeatures(
        ArchRProj = proj,
        groupBy = "Clusters_Seurat_IterativeLSI",
        $options.args
      )
    }

    markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
    sink(file = "marker_list.txt")
    for (cluster in markerList@listData) {
      cat(cluster\$name, "\n")
    }
    sink()

    # Draw heatmap: default to use all marker_genes
    if (!is.na(strtoi("$options.marker_genes"))) {
      markerGenes <- str_trim(str_split("$options.marker_genes", ","), side = "both")
      #markerGenes <- c($options.marker_genes)
    } else {
      markerGenes <- c()
      for (cluster in markerList@listData) {
        markerGenes <- c(markerGenes, cluster\$name)
      }
      sel <- min(length(markerGenes), strtoi("$options.marker_genes"))
      markerGenes <- [1:sel]
    }

    markerGenes <- unique(markerGenes)
    markerGenes_clean <- markerGenes

    ## below is to make sure genes to label are valid gene symbols in the dataset:
    all_id <- getGenes(proj)\$gene_id
    all_symbol <- getGenes(proj)\$symbol
    all_symbol_cleaned <- character(length(all_id))
    for (i in 1:length(all_id)) {
    	all_symbol_cleaned[i] <- str_remove(all_symbol[i], paste0("_", all_id[i]))
      markerGenes_clean <- str_remove(markerGenes_clean, paste0("_", all_id[i])) # not very efficient, but works
    }

    markerGenes2labeled <- sort(markerGenes_clean[markerGenes_clean %in% all_symbol_cleaned])
    markerGenes_raw <- sort(all_symbol[all_symbol_cleaned %in% markerGenes_clean])
    if (length(markerGenes2labeled) == 0) {
      message(markerGenes2labeled)
      stop("Invalid marker gene names!")
    }

    heatmapGS <- markerHeatmap(
      seMarker = markersGS,
      cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
      labelMarkers = markerGenes2labeled,
      transpose = TRUE
    )
    plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    # Plot marker genes on embeddings without imputation:
    for (embedding in names(proj@embeddings)) {
      p <- plotEmbedding(
        ArchRProj = proj,
        name = markerGenes_raw,
        imputeWeights = NULL,
        embedding = embedding,
        $options.args2
      )
      plotPDF(plotList = p, name = paste0("Plot-", embedding, "-Marker-Genes-WO-Imputation.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
    }

    # Plot marker genes on embeddings with imputation:
    proj2 <- addImputeWeights(proj)
    saveRDS(proj2, file = "proj_marker_gene.rds")
    for (embedding in names(proj@embeddings)) {
      p <- plotEmbedding(
        ArchRProj = proj2,
        name = markerGenes_raw,
        imputeWeights = getImputeWeights(proj2),
        embedding = embedding,
        $options.args2
      )
      plotPDF(plotList = p, name = paste0("Plot-", embedding, "-Marker-Genes-W-Imputation.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
    }

    # Plot: track plotting with ArchRBrowser
    clusters <- c("Clusters_Seurat_IterativeLSI", "Clusters_Scran_IterativeLSI", "Clusters_Seurat_Harmony", "Clusters_Scran_Harmony", "Clusters2_Seurat_IterativeLSI", "Clusters2_Scran_IterativeLSI", "Clusters2_Seurat_Harmony", "Clusters2_Scran_Harmony")
    for (cluster in clusters) {
      tryCatch({
        p <- plotBrowserTrack(
          ArchRProj = proj2,
          geneSymbol = markerGenes_raw,
          groupBy = cluster,
          $options.args3
        )
        plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
      },
        error=function(e) {
          message(paste0("Skipping track plotting for ", cluster, "!"))
        }
      )
    }

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
