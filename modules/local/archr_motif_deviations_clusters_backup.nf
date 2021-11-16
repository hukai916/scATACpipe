// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_MOTIF_DEVIATIONS_CLUSTERS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_motif_deviations_clusters', publish_id:'') }

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
    path "archr_project.rds", emit: archr_project

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    # Plot deviations
    proj2 <- addMotifAnnotations(ArchRProj = proj, name = "Motif", $options.args)
    proj2 <- addBgdPeaks(proj2)
    proj2 <- addDeviationsMatrix(
      ArchRProj = proj2,
      peakAnnotation = "Motif",
      force = TRUE
    )

    plotVarDev <- getVarDeviations(proj2, name = "MotifMatrix", plot = TRUE)
    plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    motifs <- c($options.motifs)
    markerMotifs <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")

    markerMotifs <- grep("z:", markerMotifs, value = TRUE)
    p <- plotGroups(ArchRProj = proj2,
      groupBy = "Clusters",
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

    # ArchR deviations
    ## Encode TFBS:
    proj2 <- addArchRAnnotations(ArchRProj = proj2, collection = "EncodeTFBS")
    proj2 <- addDeviationsMatrix(
      ArchRProj = proj2,
      peakAnnotation = "EncodeTFBS",
      force = TRUE
      )
    plotVarDev <- getVarDeviations(proj2, plot = TRUE, name = "EncodeTFBSMatrix")
    plotPDF(plotVarDev, name = "Variable-EncodeTFBS-Deviation-Scores", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    markerTFs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "EncodeTFBSMatrix")
    markerTFs <- sort(grep("z:", markerTFs, value = TRUE))
    TFnames <- stringr::str_split(stringr::str_split(markerTFs, pattern = "\\.", simplify=TRUE)[,2], pattern = "-", simplify = TRUE)[,1]
    markerTFs <- markerTFs[!duplicated(TFnames)]
    p <- plotEmbedding(
      ArchRProj = proj2,
      colorBy = "EncodeTFBSMatrix",
      name = markerTFs,
      embedding = "UMAP",
      imputeWeights = getImputeWeights(projHeme5)
      )

    ## Bulk ATAC-seq:
    proj2 <- addArchRAnnotations(ArchRProj = proj2, collection = "ATAC")
    proj2 <- addDeviationsMatrix(
      ArchRProj = proj2,
      peakAnnotation = "ATAC",
      force = TRUE
      )
    plotVarDev <- getVarDeviations(proj2, plot = TRUE, name = "ATACMatrix")
    plotPDF(plotVarDev, name = "Variable-ATAC-Deviation-Scores", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    saveRDS(proj2, file = "archr_project.rds")
    ' > run.R

    Rscript run.R

    """
}
