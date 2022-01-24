// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_EMBEDDING {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_embedding', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val archr_thread

    output:
    path "proj_embedding.rds", emit: archr_project
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_embedding", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    embedMethods <- c("UMAP", "TSNE")
    if ("Harmony" %in% names(proj@reducedDims)) {
      reducedDims  <- c("IterativeLSI", "Harmony")
    } else {
      reducedDims  <- c("IterativeLSI")
    }

    for (embedMethod in embedMethods) {
      for (reducedDim in reducedDims) {
        if (embedMethod == "UMAP") {
          proj2 <- addUMAP(
            ArchRProj = proj,
            reducedDims = reducedDim,
            name = paste0(embedMethod, "_", reducedDim),
            $options.args
          )
        } else if (embedMethod == "TSNE") {
          proj2 <- addTSNE(
            ArchRProj = proj,
            reducedDims = reducedDim,
            name = paste0(embedMethod, "_", reducedDim),
            $options.args
          )
        }
      }
    }

    saveRDS(proj2, file = "proj_embedding.rds")

    # Plotting UMAP_IterativeLSI
    p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_IterativeLSI")
    p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Seurat_IterativeLSI", embedding = "UMAP_IterativeLSI")
    p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Scran_IterativeLSI", embedding = "UMAP_IterativeLSI")
    plotPDF(p1, p2, p3, name = "Plot-UMAP-Sample-Clusters-ILSI.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    # Plotting TSNE_IterativeLSI
    p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNE_IterativeLSI")
    p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Seurat_IterativeLSI", embedding = "TSNE_IterativeLSI")
    p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Scran_IterativeLSI", embedding = "TSNE_IterativeLSI")
    plotPDF(p1, p2, p3, name = "Plot-TNSE-Sample-Clusters-ILSI.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    if ("Harmony" %in% names(proj@reducedDims)) {
      # Ploting UMAP_Harmony
      p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony")
      p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Seurat_Harmony", embedding = "UMAP_Harmony")
      p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Scran_Harmony", embedding = "UMAP_Harmony")
      plotPDF(p1, p2, p3, name = "Plot-UMAP-Sample-Clusters-Harmony.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

      # Plotting TSNE_Harmony
      p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNE_Harmony")
      p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Seurat_Harmony", embedding = "TSNE_Harmony")
      p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters_Scran_Harmony", embedding = "TSNE_Harmony")
      plotPDF(p1, p2, p3, name = "Plot-TSNE-Sample-Clusters-Harmony.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
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
      mkdir -p report_jpeg/archr_embedding
      cp -r Plots/jpeg report_jpeg/archr_embedding

    """
}
