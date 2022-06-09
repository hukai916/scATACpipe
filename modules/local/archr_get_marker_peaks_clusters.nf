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
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path archr_project
    val archr_thread

    output:
    path archr_project, emit: archr_project
    path "*_marker_peaks.rds", emit: marker_peaks
    path "*_group_names.txt", emit: group_names
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
    markerList <- getMarkers(markersPeaks, returnGR = TRUE, $options.getMarkers_cutoff)
    saveRDS(markersPeaks, file = paste0(cluster, "_marker_peaks.rds"))

    fileConn <- file(paste0(cluster, "_group_names.txt"))
    writeLines(unique(names(markerList)), fileConn)
    close(fileConn)

    tryCatch({ # use tryCatch in case no peak pass cutoff
      heatmapPeaks <- markerHeatmap(
        seMarker = markersPeaks,
        transpose = TRUE,
        $options.getMarkers_cutoff
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
        pma <- markerPlot(seMarker = markersPeaks, name = group, plotAs = "MA", $options.getMarkers_cutoff)
        pv <- markerPlot(seMarker = markersPeaks, name = group, plotAs = "Volcano", $options.getMarkers_cutoff)
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
    mkdir -p Plots/jpeg
    x=( \$(find ./Plots -name "*.pdf") )
    for item in \${x[@]+"\${x[@]}"}
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
    mkdir -p ./report_jpeg/archr_get_marker_peaks_clusters
    cp -r ./Plots/jpeg report_jpeg/archr_get_marker_peaks_clusters
    mkdir ./report_jpeg/archr_get_marker_peaks_clusters/pdf
    cp ./Plots/*.pdf report_jpeg/archr_get_marker_peaks_clusters/pdf/

    """
}
