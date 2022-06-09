include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_CLUSTERING {
  // Default to use Seurat clustering since Scran is slow and the two results are similar.
  // Check if Harmony avail, if so, use it, otherwise use IterativeLSI.
  // Only one cluster will be output at the end.

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_clustering', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path archr_project
    val filter_sample
    val filter_seurat_iLSI
    val filter_seurat_harmony
    val archr_thread

    output:
    path "proj_clustering.rds", emit: archr_project
    path "*.csv", emit: csv_cluster_matrix
    path "report_jpeg/archr_clustering", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(pheatmap)
    library(stringr)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    # get rid of undesired samples if supplied:
    if (!("$filter_sample" == "NA")) {
      filter_sample <- unique(str_trim(str_split("$filter_sample", ",")[[1]], side = "both"))
      idxPass <- which(!proj\$Sample %in% filter_sample)
      cellsPass <- proj\$cellNames[idxPass]
      proj <- proj[cellsPass,]
    }

    # first round of clustering: before removing undesired clusters
    if ("Harmony" %in% names(proj@reducedDims)) {
      # Clustering with Seurat using Harmony
      proj <- addClusters(
        input = proj,
        reducedDims = "Harmony",
        method = "Seurat",
        name = "Clusters",
        force = TRUE,
        $options.args
      )
    } else if ("IterativeLSI" %in% names(proj@reducedDims)) {
      # Clustering with Seurat using IterativeLSI
      proj <- addClusters(
        input = proj,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        force = TRUE,
        $options.args
      )
    }

    # get rid of undesired clusters if supplied and reclustering:
    if ("Harmony" %in% names(proj@reducedDims)) {
      if (!("$filter_seurat_harmony" == "NA")) {
        filter_clusters <- unique(str_trim(str_split("$filter_seurat_harmony", ",")[[1]], side = "both"))
        idxPass <- which(!proj\$Clusters %in% filter_clusters)
        cellsPass <- proj\$cellNames[idxPass]
        proj <- proj[cellsPass,]

        # second round of clustering:
        proj <- addClusters(
          input = proj,
          reducedDims = "Harmony",
          method = "Seurat",
          name = "Clusters",
          force = TRUE,
          $options.args
        )
      }
    } else if ("IterativeLSI" %in% names(proj@reducedDims)) {
      if (!("$filter_seurat_iLSI" == "NA")) {
        filter_clusters <- unique(str_trim(str_split("$filter_seurat_iLSI", ",")[[1]], side = "both"))
        idxPass <- which(!proj\$Clusters %in% filter_clusters)
        cellsPass <- proj\$cellNames[idxPass]
        proj <- proj[cellsPass,]

        # second round of clustering:
        proj <- addClusters(
          input = proj,
          reducedDims = "IterativeLSI",
          method = "Seurat",
          name = "Clusters",
          force = TRUE,
          $options.args
        )
      }
    }

    saveRDS(proj, file = "proj_clustering.rds")

    # Save text summary and heatmap summary
    cluster <- "Clusters"

    cmd <- paste0("cM <- confusionMatrix(paste0(proj\$", cluster, "), paste0(proj\$Sample))")
    eval(str2lang(cmd))
    write.csv(cM, file=paste0(cluster, "_matrix.csv"))
    cM <- cM / Matrix::rowSums(cM)
    cellheight <- 18 # 0.25 inch
    p1 <- pheatmap::pheatmap(
      mat = as.matrix(cM),
      color = paletteContinuous("whiteBlue"),
      display_numbers = TRUE,
      number_format = "%.2f",
      number_color = "red",
      fontsize_number = 8,
      border_color = "black",
      cellheight = cellheight
    )

    height <- nrow(cM) * cellheight * 1/72 + 4
    height <- min(11, height)
    # 1/72: inches per point, 11.5 inches per page; 8.5 width per page.

    plotPDF(p1, name = paste0(cluster, "_heatmap.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 7, height = height)

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
    mkdir -p ./report_jpeg/archr_clustering
    cp -r ./Plots/jpeg report_jpeg/archr_clustering
    mkdir ./report_jpeg/archr_clustering/pdf
    cp ./Plots/*.pdf report_jpeg/archr_clustering/pdf/

    """
}
