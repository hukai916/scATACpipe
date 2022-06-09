// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_EMBEDDING {
  // 4 embeddings for a single clustering ("Clusters"): UMAP/TSNE vs ILSI/Harmony
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_embedding', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path archr_project
    val archr_thread

    output:
    path "proj_embedding.rds", emit: archr_project
    path "report_jpeg/archr_embedding", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj2 <- readRDS("$archr_project", refhook = NULL)

    embedMethods <- c("UMAP", "TSNE") # TSNA may not work for some cases, use tryCatch below.
    if ("Harmony" %in% names(proj2@reducedDims)) {
      reducedDims  <- c("IterativeLSI", "Harmony")
    } else {
      reducedDims  <- c("IterativeLSI")
    }

    for (embedMethod in embedMethods) {
      for (reducedDim in reducedDims) {
        tryCatch({
          if (embedMethod == "UMAP") {
            proj2 <- addUMAP(
              ArchRProj = proj2,
              reducedDims = reducedDim,
              name = paste0(embedMethod, "_", reducedDim),
              force = TRUE,
              $options.args
            )
          } else if (embedMethod == "TSNE") {
            proj2 <- addTSNE(
              ArchRProj = proj2,
              reducedDims = reducedDim,
              name = paste0(embedMethod, "_", reducedDim),
              force = TRUE,
              $options.args
            )
          }
        },
          error=function(e) {
            message(paste0(embedMethod, "_", reducedDim, ": failed to embed!"))
          }
        )
      }
    }

    saveRDS(proj2, file = "proj_embedding.rds")

    # Plotting UMAP_IterativeLSI
    tryCatch({
      p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_IterativeLSI")
      p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_IterativeLSI")
      plotPDF(p1, p2, name = "Plot-UMAP-Sample-Clusters-ILSI.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
    },
      error=function(e) {
        message("Plotting failed: UMAP_IterativeLSI!")
      }
    )

    # Plotting TSNE_IterativeLSI
    tryCatch({
      p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNE_IterativeLSI")
      p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE_IterativeLSI")
      plotPDF(p1, p2, name = "Plot-TNSE-Sample-Clusters-ILSI.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
    },
      error=function(e) {
        message("Plotting failed: TSNE_IterativeLSI!")
      }
    )

    if ("Harmony" %in% names(proj2@reducedDims)) {
      # Ploting UMAP_Harmony
      tryCatch({
        p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony")
        p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony")
        plotPDF(p1, p2, name = "Plot-UMAP-Sample-Clusters-Harmony.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
      },
        error=function(e) {
          message("Plotting failed: UMAP_Harmony!")
        }
      )

      # Plotting TSNE_Harmony
      tryCatch({
        p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNE_Harmony")
        p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE_Harmony")
        plotPDF(p1, p2, name = "Plot-TSNE-Sample-Clusters-Harmony.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
      },
        error=function(e) {
          message("Plotting failed: TSNE_Harmony!")
        }
      )
    }

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
    mkdir -p ./report_jpeg/archr_embedding
    cp -r ./Plots/jpeg report_jpeg/archr_embedding
    mkdir ./report_jpeg/archr_embedding/pdf
    cp ./Plots/*.pdf report_jpeg/archr_embedding/pdf/

    """
}
