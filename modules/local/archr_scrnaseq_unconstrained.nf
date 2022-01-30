// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_SCRNASEQ_UNCONSTRAINED {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_scrnaseq_unconstrained', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    path obj_scrnaseq
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

    addArchRThreads(threads = $archr_thread)

    seRNA <- readRDS("$obj_scrnaseq")
    proj <- readRDS("$archr_project", refhook = NULL)

    # Perform unconstrained integration
    proj2 <- addGeneIntegrationMatrix(
      ArchRProj = proj,
      seRNA = seRNA,
      nameCell = "predictedCell_Un",
      nameGroup = "predictedGroup_Un",
      nameScore = "predictedScore_Un",
      $options.args
    )

    # Add impute weights
    proj2 <- addImputeWeights(proj2)

    # UMAP plots overlayed with gene expression values from GeneIntegrationMatrix
    markerGenes <- c($options.marker_genes)
    p1 <- plotEmbedding(
      ArchRProj = proj2,
      colorBy = "GeneIntegrationMatrix",
      name = markerGenes,
      continuousSet = "horizonExtra",
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj2)
    )

    plotPDF(plotList = p1,
      name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf",
      ArchRProj = NULL,
      addDOC = FALSE, width = 5, height = 5
    )

    # Label scATAC-seq clusters with scRNA-seq information
    cM <- confusionMatrix(proj2\$Clusters, proj2\$predictedGroup_Un)
    labelOld <- rownames(cM)
    labelNew <- colnames(cM)[apply(cM, 1, which.max)]
    proj2\$Clusters2 <- mapLabels(proj2\$Clusters, newLabels = labelNew, oldLabels = labelOld)
    p1 <- plotEmbedding(proj2, colorBy = "cellColData", name = "Clusters2")
    plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

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
