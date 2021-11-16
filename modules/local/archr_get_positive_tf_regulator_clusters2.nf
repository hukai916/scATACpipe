// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2 {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_positive_tf_regulator_clusters2', publish_id:'') }

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
    path "Plots/Plot-Tracks-With-Features.pdf", emit: plot_tracks_with_features
    path "Plots/Plot-Tracks-With-Features-Clusters2.pdf", emit: plot_tracks_with_features_clusters
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_get_positive_tf_regulator_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Clusters2")
    seZ <- seGroupMotif[rowData(seGroupMotif)\$seqnames=="z",]
    rowData(seZ)\$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
        rowMaxs(assay(seZ) - assay(seZ)[,x])
      }) %>% Reduce("cbind", .) %>% rowMaxs

    corGSM_MM <- correlateMatrices(
      ArchRProj = proj,
      useMatrix1 = "GeneScoreMatrix",
      useMatrix2 = "MotifMatrix",
      reducedDims = "IterativeLSI"
      )
    corGIM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
      )

    corGSM_MM\$maxDelta <- rowData(seZ)[match(corGSM_MM\$MotifMatrix_name, rowData(seZ)\$name), "maxDelta"]
    corGIM_MM\$maxDelta <- rowData(seZ)[match(corGIM_MM\$MotifMatrix_name, rowData(seZ)\$name), "maxDelta"]

    corGSM_MM <- corGSM_MM[order(abs(corGSM_MM\$cor), decreasing = TRUE), ]
    corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
    corGSM_MM\$TFRegulator <- "NO"
    corGSM_MM\$TFRegulator[which(corGSM_MM\$cor > 0.5 & corGSM_MM\$padj < 0.01 & corGSM_MM\$maxDelta > quantile(corGSM_MM\$maxDelta, 0.75))] <- "YES"
    sort(corGSM_MM[corGSM_MM\$TFRegulator=="YES",1])

    corGIM_MM <- corGIM_MM[order(abs(corGIM_MM\$cor), decreasing = TRUE), ]
    corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
    corGIM_MM\$TFRegulator <- "NO"
    corGIM_MM\$TFRegulator[which(corGIM_MM\$cor > 0.5 & corGIM_MM\$padj < 0.01 & corGIM_MM\$maxDelta > quantile(corGIM_MM\$maxDelta, 0.75))] <- "YES"
    sort(corGIM_MM[corGIM_MM\$TFRegulator=="YES",1])

    p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
      geom_point() +
      theme_ArchR() +
      geom_vline(xintercept = 0, lty = "dashed") +
      scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
      xlab("Correlation To Gene Score") +
      ylab("Max TF Motif Delta") +
      scale_y_continuous(
        expand = c(0,0),
        limits = c(0, max(corGSM_MM\$maxDelta)*1.05)
      )
    plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

    p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
      geom_point() +
      theme_ArchR() +
      geom_vline(xintercept = 0, lty = "dashed") +
      scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
      xlab("Correlation To Gene Expression") +
      ylab("Max TF Motif Delta") +
      scale_y_continuous(
        expand = c(0,0),
        limits = c(0, max(corGIM_MM\$maxDelta)*1.05)
      )
    plotPDF(p, name = "Plot-Tracks-With-Features-Clusters2", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

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
    mkdir -p report_jpeg/archr_get_positive_tf_regulator_clusters2
    cp -r Plots/jpeg report_jpeg/archr_get_positive_tf_regulator_clusters2

    """
}
