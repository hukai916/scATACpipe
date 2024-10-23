// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_MOTIF_ENRICHMENT_CLUSTERS2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_motif_enrichment_clusters2', publish_id:'') }
    // container "hukai916/r_archr:0.1" // must use Haibo ArchR, otherwise addMotifAnnotations problematic
    container "hukai916/r_env_v4.4.1_amd64:0.1"

    input:
    path archr_project
    path marker_test
    path markers_peaks
    path test_group
    path user_rlib
    val custom_peaks
    val species_latin_name
    val archr_thread

    output:
    path "archr_motif_enrichment_project.rds", emit: archr_project
    path "report_jpeg/archr_motif_enrichment_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    .libPaths("user_rlib") # for user installed packages

    options(timeout=10000)
    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)
    markerTest <- readRDS("$marker_test")
    markersPeaks <- readRDS("$markers_peaks")

    # Read in test group info: second line is background
    conn <- file("$test_group", open = "r")
    lines <- readLines(conn)
    useGroups <- trimws(lines[1], "both")
    bgdGroups <- trimws(lines[2], "both")
    close(conn)

    if ("$species_latin_name" == "NA") {
      proj <- addMotifAnnotations(ArchRProj = proj, name = "Motif", species = NULL, $options.args)
    } else {
      proj <- addMotifAnnotations(ArchRProj = proj, name = "Motif", species = "$species_latin_name", $options.args)
    }

    # Motif enrichment in Differential peaks:
    motifsUp <- peakAnnoEnrichment(
      seMarker = markerTest,
      ArchRProj = proj,
      peakAnnotation = "Motif",
      $options.cutoff
    )
    motifsDo <- peakAnnoEnrichment(
      seMarker = markerTest,
      ArchRProj = proj,
      peakAnnotation = "Motif",
      $options.cutoff
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

    plotPDF(ggUp, ggDo, name = paste0(useGroups, "-vs-", bgdGroups, "-Markers-Motifs-Enriched"), width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    # Motif enrichment in Marker peaks:
    enrichMotifs <- peakAnnoEnrichment(
      seMarker = markersPeaks,
      ArchRProj = proj,
      peakAnnotation = "Motif",
      $options.cutoff
    )

    tryCatch({ # use tryCatch in case no result passes cutoff
      heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
      plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)
    },
      error=function(e) {
        message(paste0("Skipping plotting motifs-enriched-marker heatmaps!"))
      }
    )

    # ArchR enrichments # is problematic when downloading using inside docker, skip it for now

    # Custom enrichment if supplied
    customPeaks <- c($custom_peaks)
    # if (!("$custom_peaks" == '') { # wont work if contain quotes
    # if (customPeaks[1] != "") { # if custome_peaks == '', NULL will be passed
    if (!(is.null(customPeaks[1]))) {
      customPeaks <- customPeaks
      proj <- addPeakAnnotations(ArchRProj = proj, regions = customPeaks, name = "Custom")
      enrichRegions <- peakAnnoEnrichment(
        seMarker = markersPeaks,
        ArchRProj = proj,
        peakAnnotation = "Custom",
        $options.cutoff
      )
      tryCatch({ # use tryCatch in case no result passes cutoff
        heatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 7, transpose = TRUE)
        plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)
      },
        error=function(e) {
          message(paste0("Skipping plotting motifs-enriched-marker heatmaps!"))
        }
      )
    }

    saveRDS(proj, file = "archr_motif_enrichment_project.rds")

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
    mkdir -p ./report_jpeg/archr_motif_enrichment_clusters2
    cp -r ./Plots/jpeg report_jpeg/archr_motif_enrichment_clusters2 || :
    mkdir ./report_jpeg/archr_motif_enrichment_clusters2/pdf
    cp ./Plots/*.pdf report_jpeg/archr_motif_enrichment_clusters2/pdf/ || :


    """
}
