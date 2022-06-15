// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_MARKER_GENE_CLUSTERS2 {
  // Find marker genes for "Clusters" only.

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_marker_gene_clusters2', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path archr_project
    val archr_thread

    output:
    path "proj_marker_gene.rds", emit: archr_project
    path "marker_list.txt", emit: marker_list
    path "Clusters2_markerList.rds", emit: markerList
    path "report_jpeg/archr_marker_gene_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(stringr)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    markersGS <- getMarkerFeatures(
      ArchRProj = proj,
      groupBy = "Clusters2",
      $options.args
    )

    markerList <- getMarkers(markersGS, $options.getMarkers_cutoff)
    sink(file = "marker_list.txt")
    for (cluster in markerList@listData) {
      cat(cluster\$name, "\n")
    }
    sink()
    saveRDS(markerList, file = paste0("Clusters2", "_markerList.rds"))

    # Draw heatmap: default to use first 3 marker_genes
    if (!("$options.marker_genes" == "default")) {
      markerGenes <- str_trim(str_split("$options.marker_genes", ",")[[1]], side = "both")
      markerGenes <- unique(markerGenes)
    } else {
      markerGenes <- c()
      for (cluster in markerList@listData) {
        markerGenes <- c(markerGenes, cluster\$name)
      }
      markerGenes <- unique(markerGenes)
      sel <- min(length(markerGenes), 3)
      markerGenes <- markerGenes[1:sel]
    }

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

    # markerGenes_raw matches with markerGenes2labelded except that the symbols match that from the dataset
    all_symbol_cleaned_unique <- unique(all_symbol_cleaned)
    all_symbol_unique <- all_symbol[match(all_symbol_cleaned_unique, all_symbol_cleaned)]
    markerGenes_raw <- sort(all_symbol_unique[all_symbol_cleaned_unique %in% markerGenes_clean])

    proj2 <- addImputeWeights(proj)
    saveRDS(proj2, file = "proj_marker_gene.rds")

    if (length(markerGenes2labeled) == 0) {
      message(markerGenes2labeled)
      message("Invalid marker gene names!")
      message("Skipping plotting!")
    } else {
      _tem <- unique(make.names(markersGS@elementMetadata@listData\$name, unique = TRUE)) # otherwise markerHeatmap won't plot label
      markersGS@elementMetadata@listData\$name <- _tem

      heatmapGS <- markerHeatmap(
        seMarker = markersGS,
        labelMarkers = markerGenes2labeled,
        binaryClusterRows = TRUE,
        clusterCols = TRUE,
        transpose = TRUE,
        $options.getMarkers_cutoff
      )
      plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 8, ArchRProj = NULL, addDOC = FALSE)

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
      # clusters <- c("Clusters_Seurat_IterativeLSI", "Clusters_Scran_IterativeLSI", "Clusters_Seurat_Harmony", "Clusters_Scran_Harmony", "Clusters2_Seurat_IterativeLSI", "Clusters2_Scran_IterativeLSI", "Clusters2_Seurat_Harmony", "Clusters2_Scran_Harmony")
      clusters <- c("Clusters2")

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
    }


    ' > run.R

    Rscript run.R

    # Convert to jpeg:
    mkdir -p Plots/jpeg
    x=( \$(find ./Plots -name "*.pdf") )
    for item in \${x[@]+"\${x[@]}"} # https://stackoverflow.com/questions/7577052/bash-empty-array-expansion-with-set-u
    do
      {
        filename=\$(basename -- "\$item")
        filename="\${filename%.*}"
        pdftoppm -jpeg -r 300 \$item ./Plots/jpeg/\$filename
        convert -append ./Plots/jpeg/\${filename}* ./Plots/jpeg/\${filename}.jpg
        rm ./Plots/jpeg/\${filename}-*.jpg
      } || {
        echo "Pdf to jpeg failed!" > bash.log
      }
    done

    # For reporting:
    mkdir -p ./report_jpeg/archr_marker_gene_clusters2
    cp -r ./Plots/jpeg report_jpeg/archr_marker_gene_clusters2
    mkdir ./report_jpeg/archr_marker_gene_clusters2/pdf
    cp ./Plots/*.pdf report_jpeg/archr_marker_gene_clusters2/pdf/

    """
}
