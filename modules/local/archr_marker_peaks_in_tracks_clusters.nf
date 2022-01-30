// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS {
  // Defaul to plot the first 10 marker genes and the first Cluster
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_marker_peaks_in_tracks_clusters', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    path marker_peaks
    path markerList
    val archr_thread

    output:
    path "Plots/Plot-Tracks-With-Features.pdf", emit: archr_tracks_with_features
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_marker_peaks_in_tracks_clusters", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(stringr)

    addArchRThreads(threads = $archr_thread)

    markersPeaks <- readRDS("$marker_peaks")
    markerList   <- readRDS("$markerList")
    proj <- readRDS("$archr_project")

    # Below is to make sure geneSymbol is subset of getGenes(proj)\$symbol
    # Draw heatmap: default to use first 10 marker_genes
    if (!("$options.marker_genes" == "default")) {
      markerGenes <- str_trim(str_split("$options.marker_genes", ",")[[1]], side = "both")
    } else {
      markerGenes <- c()
      for (cluster in markerList@listData) {
        markerGenes <- c(markerGenes, cluster\$name)
      }
      sel <- min(length(markerGenes), 10)
      markerGenes <- markerGenes[1:sel]
    }

    markerGenes <- unique(markerGenes)
    markerGenes_clean <- markerGenes

    all_id <- getGenes(proj)\$gene_id
    all_symbol <- getGenes(proj)\$symbol
    all_symbol_cleaned <- character(length(all_id))
    for (i in 1:length(all_id)) {
      all_symbol_cleaned[i] <- str_remove(all_symbol[i], paste0("_", all_id[i]))
      markerGenes_clean <- str_remove(markerGenes_clean, paste0("_", all_id[i])) # not very efficient, but works
    }

    markerGenes2labeled <- sort(markerGenes_clean[markerGenes_clean %in% all_symbol_cleaned])
    markerGenes_raw <- sort(all_symbol[all_symbol_cleaned %in% markerGenes_clean])

    all_symbol_cleaned_unique <- unique(all_symbol_cleaned)
    all_symbol_unique <- all_symbol[match(all_symbol_cleaned_unique, all_symbol_cleaned)]
    markerGenes_raw <- sort(all_symbol_unique[all_symbol_cleaned_unique %in% markerGenes_clean])

    if (length(markerGenes2labeled) == 0) {
      message(markerGenes2labeled)
      message("Invalid marker gene names!")
      message("Skipping plotting tracks-with-features!")
    } else {
      # Below is to decide cluster_name, default to use the first one
      if ("$options.cluster_name" == "default") {
        cluster_name <- 1
      } else {
        cluster_name <- "$options.cluster_name"
      }

      p <- plotBrowserTrack(
        ArchRProj = proj,
        groupBy = "Clusters",
        geneSymbol = markerGenes_raw,
        features =  getMarkers(markersPeaks, cutOff = "$options.getMarkers_cutoff", returnGR = TRUE)[cluster_name],
        $options.args
      ) # if p == 0, the pdf will be empty, and the converting to jpeg is problematic.

      plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)
    }


    ' > run.R

    Rscript run.R

    # Convert to jpeg:
    mkdir -p Plots/jpeg
    x=( \$(find ./Plots -name "*.pdf") )
    for item in "\${x[@]}"
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
    mkdir -p report_jpeg/archr_marker_peaks_in_tracks_clusters
    cp -r Plots/jpeg report_jpeg/archr_marker_peaks_in_tracks_clusters

    """
}
