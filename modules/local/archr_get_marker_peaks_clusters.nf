// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_GET_MARKER_PEAKS_CLUSTERS {
  //
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_marker_peaks_clusters', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val archr_thread

    output:
    path "Clusters_marker_peaks.rds", emit: archr_marker_peaks
    path "*_group_names.txt", emit: group_names
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_get_marker_peaks_clusters", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    cluster <- "Clusters"
    markersPeaks <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = cluster,
      $options.args
    )
    saveRDS(markersPeaks, file = paste0(cluster, "_marker_peaks.rds"))

    markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff", returnGR = TRUE)
    fileConn <- file(paste0(cluster, "_group_names.txt"))
    writeLines(unique(names(markerList)), fileConn)
    close(fileConn)

    tryCatch({ # use tryCatch in case no peak pass cutoff
      heatmapPeaks <- markerHeatmap(
        seMarker = markersPeaks,
        cutOff = "$options.cutoff",
        transpose = TRUE
      )
      # height <- nrow(cM) * cellheight * 1/72 + 4
      # height <- min(11, height)
      plotPDF(heatmapPeaks, name = paste0(cluster, "-Peak-Marker-Heatmap"), width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)
    },
      error=function(e) {
        message(paste0("Skipping plotting heatmaps!"))
      }
    )

    # Plot MA and Vocalno plots for each group:
    for (group in unique(names(markerList))) {
      tryCatch({
        pma <- markerPlot(seMarker = markersPeaks, cutOff = "$options.cutoff", name = group, plotAs = "MA")
        pv <- markerPlot(seMarker = markersPeaks, cutOff = "$options.cutoff", name = group, plotAs = "Volcano")
        plotPDF(pma, pv, name = paste0(cluster, "-", group, "-Markers-MA-Volcano"), width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)
      },
        error=function(e) {
          message(paste0("Skipping plotting Vocalno plots for ", group, "!"))
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
    mkdir -p report_jpeg/archr_get_marker_peaks_clusters
    cp -r Plots/jpeg report_jpeg/archr_get_marker_peaks_clusters

    """
}
