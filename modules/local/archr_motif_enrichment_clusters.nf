// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_MOTIF_ENRICHMENT_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_motif_enrichment_clusters', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    path marker_test
    path markers_peaks
    val useGroups
    val bgdGroups
    val custom_peaks
    val archr_thread

    output:
    path "archr_motif_enrichment_project.rds", emit: archr_project
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_motif_enrichment_clusters", emit: report

    script:

    """
    echo '
    library(ArchR)

    options(timeout=10000)
    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)
    markerTest <- readRDS("$marker_test")
    markersPeaks <- readRDS("$markers_peaks")

    proj2 <- addMotifAnnotations(ArchRProj = proj, name = "Motif", $options.args)

    # Motif enrichment in Differential peaks:
    motifsUp <- peakAnnoEnrichment(
      seMarker = markerTest,
      ArchRProj = proj2,
      peakAnnotation = "Motif",
      cutOff = "$options.cutoff"
    )
    motifsDo <- peakAnnoEnrichment(
      seMarker = markerTest,
      ArchRProj = proj2,
      peakAnnotation = "Motif",
      cutOff = "$options.cutoff"
    )
    df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df <- df[order(df\$mlog10Padj, decreasing = TRUE),]
    df\$rank <- seq_len(nrow(df))
    ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
            data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
            size = 1.5,
            nudge_x = 2,
            color = "black"
      ) + theme_ArchR() +
      ylab("-log10(P-adj) Motif Enrichment") +
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))

    df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
    df <- df[order(df\$mlog10Padj, decreasing = TRUE),]
    df\$rank <- seq_len(nrow(df))
    ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
            data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
            size = 1.5,
            nudge_x = 2,
            color = "black"
      ) + theme_ArchR() +
      ylab("-log10(FDR) Motif Enrichment") +
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))

    plotPDF(ggUp, ggDo, name = paste0("$useGroups", "-vs-", "$bgdGroups", "-Markers-Motifs-Enriched"), width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    # Motif enrichment in Marker peaks:
    enrichMotifs <- peakAnnoEnrichment(
      seMarker = markersPeaks,
      ArchRProj = proj2,
      peakAnnotation = "Motif",
      cutOff = "$options.cutoff"
    )
    heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
    plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    # ArchR enrichments # is problematic when downloading using inside docker, skip it for now

    # Custom enrichment if supplied
    customPeaks <- c($custom_peaks)
    # if (!("$custom_peaks" == '') { # wont work if contain quotes
    # if (customPeaks[1] != "") { # if custome_peaks == '', NULL will be passed
    if (!(is.null(customPeaks[1]))) {
      customPeaks <- customPeaks
      proj2 <- addPeakAnnotations(ArchRProj = proj2, regions = customPeaks, name = "Custom")
      enrichRegions <- peakAnnoEnrichment(
        seMarker = markersPeaks,
        ArchRProj = proj2,
        peakAnnotation = "Custom",
        cutOff = "$options.cutoff"
      )
      heatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 7, transpose = TRUE)
      plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)
    }

    saveRDS(proj2, file = "archr_motif_enrichment_project.rds")

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
    mkdir -p report_jpeg/archr_motif_enrichment_clusters
    cp -r Plots/jpeg report_jpeg/archr_motif_enrichment_clusters

    """
}
