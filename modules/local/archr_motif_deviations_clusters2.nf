// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_MOTIF_DEVIATIONS_CLUSTERS2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_motif_deviations_clusters2', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val custom_peaks
    val archr_thread

    output:
    path "archr_motif_deviation_project.rds", emit: archr_project
    path "motif_names.txt", emit: motif_names
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_motif_deviations_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    # Plot deviations
    # proj2 <- addMotifAnnotations(ArchRProj = proj, name = "Motif", $options.args)
    proj2 <- addBgdPeaks(proj, force = TRUE)
    proj2 <- addDeviationsMatrix(
      ArchRProj = proj2,
      peakAnnotation = "Motif",
      force = TRUE
    )

    plotVarDev <- getVarDeviations(proj2, name = "MotifMatrix", plot = TRUE)
    VarDev     <- getVarDeviations(proj2, name = "MotifMatrix", plot = FALSE)
    sink(file = "motif_names.txt")
    for (motif in VarDev\$name) {
      cat(motif, "\n")
    }
    sink()

    plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    motifs <- c($options.motifs)
    markerMotifs <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")

    markerMotifs <- grep("z:", markerMotifs, value = TRUE)
    p <- plotGroups(ArchRProj = proj2,
      groupBy = "Clusters2",
      colorBy = "MotifMatrix",
      name = markerMotifs,
      imputeWeights = getImputeWeights(proj2)
      )
    plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    # With UMAP embedding:
    p <- plotEmbedding(
      ArchRProj = proj2,
      colorBy = "MotifMatrix",
      name = sort(markerMotifs),
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj2)
      )
    plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-UMAP-Embedding", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    # With gene scores overlay on the UMAP
    markerRNA <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
    p <- plotEmbedding(
      ArchRProj = proj2,
      colorBy = "GeneScoreMatrix",
      name = sort(markerRNA),
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj2)
      )
    plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-UMAP-Embedding-w-Gene-Scores", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    # With gene expression overlay on the UMAP: for clusters2 only
    # ArchR deviations, skipped before fixing the Nextflow download.file problem.
    markerRNA <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
    p <- plotEmbedding(
      ArchRProj = proj2,
      colorBy = "GeneIntegrationMatrix",
      name = sort(markerRNA),
      embedding = "UMAP",
      continuousSet = "blueYellow",
      imputeWeights = getImputeWeights(proj2)
    )
    plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-UMAP-Embedding-w-Gene-Expression", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    # Custom enrichment if supplied
    customPeaks <- c($custom_peaks)
    # if (!("$custom_peaks" == '') { # wont work if contain quotes
    # if (customPeaks[1] != "") {
    if (!(is.null(customPeaks[1]))) {
      customPeaks <- customPeaks
      # proj2 <- addPeakAnnotations(ArchRProj = proj2, regions = customPeaks, name = "Custom")

      proj2 <- addDeviationsMatrix(
        ArchRProj = proj2,
        peakAnnotation = "Custom",
        force = TRUE
        )
      plotVarDev <- getVarDeviations(proj2, plot = TRUE, name = "CustomMatrix")
      plotPDF(plotVarDev, name = "Variable-Custom-Deviation-Scores", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

      markerCustom <- getFeatures(proj2, useMatrix = "CustomMatrix")
      markerCustom <- sort(grep("z:", markerCustom, value = TRUE))

      p <- plotEmbedding(
        ArchRProj = proj2,
        colorBy = "CustomMatrix",
        name = markerCustom,
        embedding = "UMAP",
        imputeWeights = getImputeWeights(proj2)
        )

      plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-UMAP-Embedding-w-z-scores", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)
    }

    saveRDS(proj2, file = "archr_motif_deviation_project.rds")
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
    mkdir -p report_jpeg/archr_motif_deviations_clusters2
    cp -r Plots/jpeg report_jpeg/archr_motif_deviations_clusters2

    """
}
