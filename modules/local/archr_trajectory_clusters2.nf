// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_TRAJECTORY_CLUSTERS2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_trajectory_clusters2', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val archr_thread

    output:
    path "archr_trajectory_project.rds", emit: archr_project
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "archr_trajectory_clusters2", emit: res_dir
    path "report_jpeg/archr_trajectory_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    if ("$options.trajectory_groups" == "default") {
      # Use the first 3 groups, no biological meaning
      clusters2 <- unique(proj\$Clusters2)
      trajectory <- clusters2[1:min(length(clusters2), 3)]
    } else {
      trajectory <- str_trim(str_split("$options.trajectory_groups", ",")[[1]], side = "both")
    }

    if (!(all(trajectory %in% clusters2))) {
      stop("Not all trajectory_groups are valid clusters2 group names!")
    }

    add_trajectory <- function(embedding) {
      trajectory_name <- paste0("$options.trajectory_name", "-", embedding)
      proj <- addTrajectory(
                  ArchRProj = proj,
                  name = trajectory_name,
                  groupBy = "Clusters2",
                  trajectory = trajectory,
                  embedding = embedding,
                  force = TRUE
               )

      p <- plotTrajectory(proj, embedding = embedding, trajectory = trajectory_name, colorBy = "cellColData", name = trajectory_name)
      plotPDF(p, name = paste0("Plot-Traj-", embedding, ".pdd"), ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

      # overlay the specified gene onto the embedding
      p1 <- plotTrajectory(proj, embedding = embedding, trajectory = trajectory_name, colorBy = "GeneScoreMatrix", name = trajectory_name, continuousSet = "horizonExtra")
      p2 <- plotTrajectory(proj, embedding = embedding, trajectory = trajectory_name, colorBy = "GeneIntegrationMatrix", name = trajectory_name, continuousSet = "blueYellow")

      plotPDF(p1, name = paste0("Plot-Traj-", embedding, "-w-GeneScore.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)
      plotPDF(p2, name = paste0("Plot-Traj-", embedding, "-w-GeneExpression.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

      # Plot pseudo-time heatmaps for motifs, gene scores, gene expression and peak accessibility
      trajMM  <- getTrajectory(ArchRProj = proj, name = trajectory_name, useMatrix = "MotifMatrix", log2Norm = FALSE)
      p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

      trajGSM <- getTrajectory(ArchRProj = proj, name = trajectory_name, useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
      p2 <- trajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"))

      trajGIM <- getTrajectory(ArchRProj = proj, name = trajectory_name, useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
      p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))

      trajPM  <- getTrajectory(ArchRProj = proj, name = trajectory_name, useMatrix = "PeakMatrix", log2Norm = TRUE)
      p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))

      plotPDF(p1, p2, p3, p4, name = paste0("Plot-Traj-", embedding, "-Heatmaps.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 6, height = 8)

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
      plotPDF(ht1 + ht2, name = paste0("Plot-Traj-", embedding, "-Paired-Heatmaps-w-GeneScore.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 6, height = 8)

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
      plotPDF(ht1 + ht2, name = paste0("Plot-Traj-", embedding, "Paired-Heatmaps-w-GeneExpression.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 6, height = 8)

      return(proj)
    }

    # Plot overlay to every embedding:
    for (embedding in names(proj@embeddings)) {
      proj <- add_trajectory(embedding)
    }

    saveRDS(proj, file = "archr_trajectory_project.rds")

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
    mkdir -p report_jpeg/archr_trajectory_clusters2
    cp -r Plots/jpeg report_jpeg/archr_trajectory_clusters2

    """
}
