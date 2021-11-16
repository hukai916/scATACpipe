// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_GET_MARKER_PEAKS_CLUSTERS2 {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_marker_peaks_clusters2', publish_id:'') }

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
    path "marker_peaks.rds", emit: archr_marker_peaks
    path "Plots/Peak-Marker-Heatmap.pdf", emit: archr_peak_marker_heatmap
    path "group_names.txt", emit: group_names
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_get_marker_peaks_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    markersPeaks <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = "Clusters2",
      $options.args
    )

    saveRDS(markersPeaks, file = "marker_peaks.rds")

    markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff", returnGR = TRUE)
    fileConn <- file("group_names.txt")
    writeLines(unique(names(markerList)), fileConn)
    close(fileConn)

    heatmapPeaks <- markerHeatmap(
      seMarker = markersPeaks,
      cutOff = "$options.cutoff",
      transpose = TRUE
    )

    # height <- nrow(cM) * cellheight * 1/72 + 4
    # height <- min(11, height)
    plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    # Plot MA and Vocalno plots for each group:
    for (group in unique(names(markerList))) {
      pma <- markerPlot(seMarker = markersPeaks, cutOff = "$options.cutoff", name = group, plotAs = "MA")
      pv <- markerPlot(seMarker = markersPeaks, cutOff = "$options.cutoff", name = group, plotAs = "Volcano")
      plotPDF(pma, pv, name = paste0(group, "-Markers-MA-Volcano"), width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)
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
    mkdir -p report_jpeg/archr_get_marker_peaks_clusters2
    cp -r Plots/jpeg report_jpeg/archr_get_marker_peaks_clusters2

    """
}
