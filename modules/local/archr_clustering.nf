// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_CLUSTERING {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_clustering', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val filter_seurat_iLSI
    val filter_seurat_harmony
    val filter_scran_iLSI
    val filter_scran_harmony
    val archr_thread

    output:
    path "proj_clustering.rds", emit: archr_project
    path "*.csv", emit: csv_cluster_matrix
    path "Plots/*.pdf", emit: pdf_cluster_heatmap
    path "Plots/jpeg", emit: jpeg
    path "report_jpeg/archr_clustering", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(pheatmap)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    # Clustering with Seurat using IterativeLSI
    proj2 <- addClusters(
      input = proj,
      reducedDims = "IterativeLSI",
      method = "Seurat",
      name = "Clusters_Seurat_IterativeLSI",
      $options.args
    )
    # Clustering with Scran using IterativeLSI
    proj2 <- addClusters(
      input = proj2,
      reducedDims = "IterativeLSI",
      method = "scran",
      name = "Clusters_Scran_IterativeLSI",
      $options.args2
    )

    if ("Harmony" %in% names(proj@reducedDims)) {
      # Clustering with Seurat using Harmony
      proj2 <- addClusters(
        input = proj,
        reducedDims = "Harmony",
        method = "Seurat",
        name = "Clusters_Seurat_Harmony",
        $options.args
      )
      # Clustering with Scran using Harmony
      proj2 <- addClusters(
        input = proj2,
        reducedDims = "Harmony",
        method = "scran",
        name = "Clusters_Scran_Harmony",
        $options.args2
      )
    }

    # get rid of undesired clusters if supplied:
    if (!($filter_seurat_iLSI == "NA")) {
      idxPass <- which(!proj2\$Clusters_Seurat_IterativeLSI %in% c($filter_seurat_iLSI))
      cellsPass <- proj2\$cellNames[idxPass]
      proj2 <- proj2[cellsPass,]
    }
    if (!($filter_scran_iLSI == "NA")) {
      idxPass <- which(!proj2\$Clusters_Scran_IterativeLSI %in% c($filter_scran_iLSI))
      cellsPass <- proj2\$cellNames[idxPass]
      proj2 <- proj2[cellsPass,]
    }

    if ("Harmony" %in% names(proj@reducedDims)) {
      if (!($filter_seurat_harmony == "NA")) {
        idxPass <- which(!proj2\$Clusters_Seurat_Harmony %in% c($filter_seurat_harmony))
        cellsPass <- proj2\$cellNames[idxPass]
        proj2 <- proj2[cellsPass,]
      }
      if (!($filter_scran_harmony == "NA")) {
        idxPass <- which(!proj2\$Clusters_Scran_Harmony %in% c($filter_scran_harmony))
        cellsPass <- proj2\$cellNames[idxPass]
        proj2 <- proj2[cellsPass,]
      }
    }

    saveRDS(proj2, file = "proj_clustering.rds")

    # Save text summary and heatmap summary
    if ("Harmony" %in% names(proj@reducedDims)) {
      clusters <- c("Clusters_Seurat_IterativeLSI", "Clusters_Seurat_Harmony", "Clusters_Scran_IterativeLSI", "Clusters_Scran_Harmony")
    } else {
      clusters <- c("Clusters_Seurat_IterativeLSI", "Clusters_Scran_IterativeLSI")
    }

    for (cluster in clusters) {
      cmd <- paste0("cM <- confusionMatrix(paste0(proj2\$", cluster, "), paste0(proj2\$Sample))")
      eval(str2lang(cmd))
      write.csv(cM, file=paste0(cluster, "_matrix.csv"))
      cM <- cM / Matrix::rowSums(cM)
      cellheight <- 18 # 0.25 inch
      p1 <- pheatmap::pheatmap(
        mat = as.matrix(cM),
        color = paletteContinuous("whiteBlue"),
        border_color = "black",
        cellheight = cellheight
      )
      plotPDF(p1, name = paste0(cluster, "_heatmap.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 7, height = height)
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
    mkdir -p report_jpeg/archr_clustering
    cp -r Plots/jpeg report_jpeg/archr_clustering

    """
}
