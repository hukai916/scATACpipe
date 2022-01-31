// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_SCRNASEQ_UNCONSTRAINED {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_scrnaseq_unconstrained', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    path obj_scrnaseq
    path archr
    val archr_thread

    output:
    path "proj_scrnaseq_unconstrained.rds", emit: archr_project
    path "cell_type_scRNA.txt", emit: cell_type_scrna
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_scrnaseq_unconstrained", emit: report

    script:

    """
    echo '
    library(ArchR)
    devtools::load_all("ArchR") # to use .getFeatureDF() function

    addArchRThreads(threads = $archr_thread)

    seRNA <- readRDS("$obj_scrnaseq")
    proj <- readRDS("$archr_project", refhook = NULL)

    # Perform unconstrained integration
    if ("Harmony" %in% names(proj@reducedDims)) {
      reducedDims <- "Harmony"
    } else if ("IterativeLSI" %in% names(proj@reducedDims)) {
      reducedDims <- "IterativeLSI"
    }

    proj2 <- addGeneIntegrationMatrix(
      ArchRProj = proj,
      seRNA = seRNA,
      nameCell = "predictedCell_Un",
      nameGroup = "predictedGroup_Un",
      nameScore = "predictedScore_Un",
      useMatrix = "GeneScoreMatrix",
      matrixName = "GeneIntegrationMatrix",
      reducedDims = reducedDims,
      addToArrow = TRUE,
      force = TRUE,
      $options.args
    )

    # Add impute weights
    proj2 <- addImputeWeights(proj2)

    # Embedding plots overlayed with gene expression values from GeneIntegrationMatrix
    # markerGenes: default to use first 10 marker_genes inferred from "Clusters"
    if (!("$options.marker_genes" == "default")) {
      markerGenes <- str_trim(str_split("$options.marker_genes", ",")[[1]], side = "both")
    } else {
      markersGS <- getMarkerFeatures(
        ArchRProj = proj,
        groupBy = "Clusters",
        useMatrix = "GeneScoreMatrix",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
      )
      markerList <- getMarkers(markersGS)

      markerGenes <- c()
      for (cluster in markerList@listData) {
        markerGenes <- c(markerGenes, cluster\$name)
      }
      sel <- min(length(markerGenes), 10)
      markerGenes <- markerGenes[1:sel]
    }
    # markerGenes must also be a subset of .getFeatureDF(getArrowFiles(proj2), "GeneIntegrationMatrix")\$name
    geneDF <- .getFeatureDF(getArrowFiles(proj2), "GeneScoreMatrix")
    geneDF <- geneDF[geneDF\$name %in% rownames(seRNA), , drop = FALSE]
    markerGenes <- markerGenes[markerGenes %in% geneDF\$name]

    # Plotting for embedding: can only choose one embedding, default to use UMAP.
    embedding <- paste0("UMAP_", reducedDims)

    p1 <- plotEmbedding(
      ArchRProj = proj2,
      colorBy = "GeneIntegrationMatrix",
      name = markerGenes,
      continuousSet = "horizonExtra",
      embedding = embedding,
      imputeWeights = getImputeWeights(proj2)
    )

    plotPDF(plotList = p1,
      name = paste0("Plot-", embedding, "-Marker-Genes-RNA-W-Imputation.pdf"),
      ArchRProj = NULL,
      addDOC = FALSE, width = 5, height = 5
    )

    # Label scATAC-seq clusters with scRNA-seq information
    cM <- confusionMatrix(proj2\$Clusters, proj2\$predictedGroup_Un)
    labelOld <- rownames(cM)
    labelNew <- colnames(cM)[apply(cM, 1, which.max)]
    proj2\$Clusters2 <- mapLabels(proj2\$Clusters, newLabels = labelNew, oldLabels = labelOld)
    p1 <- plotEmbedding(proj2, colorBy = "cellColData", name = "Clusters2")
    plotPDF(p1, name = paste0("Plot-", embedding, "-Remap-Clusters.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    saveRDS(proj2, file = "proj_scrnaseq_unconstrained.rds")
    fileConn <- file("cell_type_scRNA.txt")
    writeLines(unique(proj2\$predictedGroup_Un), fileConn)
    close(fileConn)

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
    mkdir -p report_jpeg/archr_scrnaseq_unconstrained
    cp -r Plots/jpeg report_jpeg/archr_scrnaseq_unconstrained

    """
}
