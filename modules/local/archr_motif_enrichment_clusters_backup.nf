// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_MOTIF_ENRICHMENT_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_motif_enrichment_clusters', publish_id:'') }

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
    path marker_test
    path markers_peaks
    val useGroups
    val bgdGroups
    val custom_peaks

    output:
    path "archr_project.rds", emit: archr_project
    path "Plots/*-vs-*-Markers-Motifs-Enriched.pdf", emit: markers_motifs_enriched
    path "Plots/Motifs-Enriched-Marker-Heatmap.pdf", emit: motifs_enriched_marker_heatmap
    path "Plots/EncodeTFBS-Enriched-Marker-Heatmap.pdf", emit: encodetfbs_marker_heatmap
    path "Plots/ATAC-Enriched-Marker-Heatmap.pdf", emit: atac_marker_heatmap
    path "Plots/Codex-Enriched-Marker-Heatmap.pdf", emit: codex_marker_heatmap

    script:

    """
    echo '
    options(timeout=10000)
    library(ArchR)

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

    # ArchR enrichments
    proj2 <- addArchRAnnotations(ArchRProj = proj2, collection = "EncodeTFBS")
    enrichEncode <- peakAnnoEnrichment(
      seMarker = markersPeaks,
      ArchRProj = proj2,
      peakAnnotation = "EncodeTFBS",
      cutOff = "$options.cutoff"
    )
    heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
    plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    proj2 <- addArchRAnnotations(ArchRProj = proj2, collection = "ATAC")
    enrichATAC <- peakAnnoEnrichment(
      seMarker = markersPeaks,
      ArchRProj = proj2,
      peakAnnotation = "ATAC",
      cutOff = "$options.cutoff"
    )
    heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
    plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    proj2 <- addArchRAnnotations(ArchRProj = proj2, collection = "Codex")
    enrichCodex <- peakAnnoEnrichment(
      seMarker = markersPeaks,
      ArchRProj = proj2,
      peakAnnotation = "Codex",
      cutOff = "$options.cutoff"
    )
    heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
    plotPDF(heatmapCodex, name = "Codex-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)

    # Custom enrichment if supplied
    customPeaks <- c($custom_peaks)
    # if (!("$custom_peaks" == '') { # wont work if contain quotes
    if (customPeaks[1] != "") {
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

    saveRDS(proj2, file = "archr_project.rds")

    ' > run.R

    Rscript run.R

    """
}
