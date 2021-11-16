// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_EMBEDDING {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_embedding', publish_id:'') }

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
    path "proj_embedding.rds", emit: archr_project
    path "Plots/Plot-UMAP-Sample-Clusters.pdf", emit: pdf_umap_sample_clusters
    path "Plots/Plot-UMAP-Sample-ScranClusters.pdf", emit: pdf_umap_sample_scranclusters
    path "Plots/Plot-tSNE-Sample-Clusters.pdf", emit: pdf_tsne_sample_clusters
    path "Plots/Plot-tSNE-Sample-ScranClusters.pdf", emit: pdf_tsne_sample_scranclusters
    path "Plots/Plot-UMAP2Harmony-Sample-Clusters.pdf", emit: pdf_umap2harmony_sample_clusters
    path "Plots/Plot-TSNE2Harmony-Sample-Clusters.pdf", emit: pdf_tsne2harmony_sample_clusters
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_embedding", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- addUMAP(
      ArchRProj = proj,
      reducedDims = "IterativeLSI",
      name = "UMAP",
      $options.args
    )

    proj2 <- addTSNE(
      ArchRProj = proj2,
      reducedDims = "IterativeLSI",
      name = "TSNE",
      $options.args2
    )

    proj2 <- addUMAP(
      ArchRProj = proj2,
      reducedDims = "Harmony",
      name = "UMAPHarmony",
      $options.args
    )

    proj2 <- addTSNE(
      ArchRProj = proj2,
      reducedDims = "Harmony",
      name = "TSNEHarmony",
      $options.args2
    )

    saveRDS(proj2, file = "proj_embedding.rds")

    # Plotting for UMAP with seurat clustering:
    p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
    p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
    plotPDF(p1, p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    # Plotting for UMAP with scran clustering:
    p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
    p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "ScranClusters", embedding = "UMAP")
    plotPDF(p3, p4, name = "Plot-UMAP-Sample-ScranClusters.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)


    # Plotting for tSNE with seurat clustering:
    p11 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
    p22 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
    plotPDF(p11, p22, name = "Plot-tSNE-Sample-Clusters.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    # Plotting for tSNE with scran clustering:
    p33 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
    p44 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "ScranClusters", embedding = "TSNE")
    plotPDF(p33, p44, name = "Plot-tSNE-Sample-ScranClusters.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)


    # Plotting for batch correctd UMAP with seurat clustering:
    p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
    p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
    plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    # Plotting for batch corrected tSNE with seurat clustering:
    p33 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
    p44 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
    plotPDF(p11,p22,p33,p44, name = "Plot-TSNE2Harmony-Sample-Clusters.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

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
