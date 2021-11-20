// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2 {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_marker_peaks_in_tracks_clusters2', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    path marker_peaks
    val gene_symbol
    val cluster_name

    output:
    path "Plots/Plot-Tracks-With-Features.pdf", emit: archr_tracks_with_features
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_marker_peaks_in_tracks_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project")
    markersPeaks <- readRDS("$marker_peaks")

    p <- plotBrowserTrack(
      ArchRProj = proj,
      groupBy = "Clusters2",
      geneSymbol = c("$gene_symbol"),
      features =  getMarkers(markersPeaks, cutOff = "$options.cutoff", returnGR = TRUE)["$cluster_name"],
      $options.args
    )

    plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

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
    mkdir -p report_jpeg/archr_marker_peaks_in_tracks_clusters2
    cp -r Plots/jpeg report_jpeg/archr_marker_peaks_in_tracks_clusters2

    """
}
