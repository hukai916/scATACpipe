// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_MARKER_PEAKS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_marker_peaks', publish_id:'') }

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
    val groupby

    output:
    path "proj_call_peaks.rds", emit: archr_project

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes().intdiv(task.cpus)}" : ''
    // Ref: https://gitter.im/nextflow-io/nextflow?at=5a4f8f01ce68c3bc7480d7c5


    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)
    if ($groupby == "Clusters") {
      # call marker peaks on Clusters only

      markerPeaks <- getMarkerFeatures(
        ArchRProj = proj,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters",
        $options.args
      )

      markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff")
      markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff", returnGR = TRUE)
      heatmapPeaks <- markerHeatmap(
        seMarker = markersPeaks,
        cutOff = "$options.cutoff",
        transpose = TRUE
      )

      # height <- nrow(cM) * cellheight * 1/72 + 4
      # height <- min(11, height)

      plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-Clusters", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    } else if ($groupby == "Clusters2") {
      # call marker peaks on both Clusters2 and Clusters
      markerPeaks <- getMarkerFeatures(
        ArchRProj = proj,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters",
        $options.args
      )

      markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff")
      markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff", returnGR = TRUE)
      heatmapPeaks <- markerHeatmap(
        seMarker = markersPeaks,
        cutOff = "$options.cutoff",
        transpose = TRUE
      )

      # height <- nrow(cM) * cellheight * 1/72 + 4
      # height <- min(11, height)

      plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-Clusters", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

      markerPeaks <- getMarkerFeatures(
        ArchRProj = proj,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters2",
        $options.args
      )

      markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff")
      markerList <- getMarkers(markersPeaks, cutOff = "$options.cutoff", returnGR = TRUE)
      heatmapPeaks <- markerHeatmap(
        seMarker = markersPeaks,
        cutOff = "$options.cutoff",
        transpose = TRUE
      )

      # height <- nrow(cM) * cellheight * 1/72 + 4
      # height <- min(11, height)

      plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-Clusters2", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)
    }

    ' > run.R

    Rscript run.R

    """
}
