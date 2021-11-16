// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_TRAJECTORY_CLUSTERS2 {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_trajectory_clusters2', publish_id:'') }

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
    val trajectory_groups

    output:
    path "archr_trajectory_project.rds", emit: archr_project
    path "Plots/Plot-Traj-UMAP.pdf", emit: plot_traj_umap
    path "Plots/Plot-Traj-Heatmaps.pdf", emit: plot_traj_heatmap
    path "Plots/Plot-Traj-Paired-Heatmaps-w-GeneScore.pdf", emit: plot_traj_paired_heatmaps_w_genescore
    path "Plots/Plot-Traj-Paired-Heatmaps-w-GeneExpression.pdf", emit: plot_traj_paired_heatmaps_w_geneexpression
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "archr_trajectory_clusters2", emit: res_dir
    path "report_jpeg/archr_trajectory_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    trajectory = c($trajectory_groups)

    proj2 <- addTrajectory(
      ArchRProj = proj,
      name = "Trajectory",
      groupBy = "Clusters2",
      trajectory = trajectory,
      embedding = "UMAP",
      force = TRUE
      )

    saveRDS(proj2, file = "archr_trajectory_project.rds")

    # Plot UMAP
    p <- plotTrajectory(proj2, trajectory = "Trajectory", colorBy = "cellColData", name = "Trajectory")
    plotPDF(p, name = "Plot-Traj-UMAP.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    if (!($options.gene_to_color == "")) {
      # overlay the specified gene onto the UMAP
      p1 <- plotTrajectory(proj2, trajectory = "Trajectory", colorBy = "GeneScoreMatrix", name = $options.gene_to_color, continuousSet = "horizonExtra")
      p2 <- plotTrajectory(proj2, trajectory = "Trajectory", colorBy = "GeneIntegrationMatrix", name = $options.gene_to_color, continuousSet = "blueYellow")

      plotPDF(p1, name = "Plot-Traj-UMAP-w-GeneScore.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
      plotPDF(p2, name = "Plot-Traj-UMAP-w-GeneExpression.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
      }

    # Plot pseudo-time heatmaps for motifs, gene scores, gene expression and peak accessibility
    trajMM  <- getTrajectory(ArchRProj = proj2, name = "Trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
    p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

    trajGSM <- getTrajectory(ArchRProj = proj2, name = "Trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
    p2 <- trajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"))

    trajGIM <- getTrajectory(ArchRProj = proj2, name = "Trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
    p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))

    trajPM  <- getTrajectory(ArchRProj = proj2, name = "Trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
    p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))

    plotPDF(p1, p2, p3, p4, name = "Plot-Traj-Heatmaps.pdf", ArchRProj = NULL, addDOC = FALSE, width = 6, height = 8)

    # Integrative pseudo-time analyses
    corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
    trajGSM2 <- trajGSM[corGSM_MM[[1]]\$name1, ]
    trajMM2 <- trajMM[corGSM_MM[[1]]\$name2, ]
    trajCombined <- trajGSM2
    assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

    combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
    rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
    ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
    ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
    plotPDF(ht1 + ht2, name = "Plot-Traj-Paired-Heatmaps-w-GeneScore.pdf", ArchRProj = NULL, addDOC = FALSE, width = 6, height = 8)

    corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
    trajGIM2 <- trajGIM[corGIM_MM[[1]]\$name1, ]
    trajMM2 <- trajMM[corGIM_MM[[1]]\$name2, ]
    trajCombined <- trajGIM2
    assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

    rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
    combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
    rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
    ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
    ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
    plotPDF(ht1 + ht2, name = "Plot-Traj-Paired-Heatmaps-w-GeneExpression.pdf", ArchRProj = NULL, addDOC = FALSE, width = 6, height = 8)

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
    mkdir -p report_jpeg/archr_trajectory_clusters2
    cp -r Plots/jpeg report_jpeg/archr_trajectory_clusters2

    """
}
