// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_positive_tf_regulator_clusters2', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path archr_project
    val archr_thread

    output:
    path "report_jpeg/archr_get_positive_tf_regulator_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(ggrepel)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    tryCatch({
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

      corGSM_MM\$maxDelta <- rowData(seZ)[match(corGSM_MM\$MotifMatrix_name, rowData(seZ)\$name), "maxDelta"]

      corGSM_MM <- corGSM_MM[order(abs(corGSM_MM\$cor), decreasing = TRUE), ]
      corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
      corGSM_MM\$TFRegulator <- "NO"
      corGSM_MM\$TFRegulator[which(corGSM_MM\$cor > 0.5 & corGSM_MM\$padj < 0.01 & corGSM_MM\$maxDelta > quantile(corGSM_MM\$maxDelta, 0.75))] <- "YES"
      sort(corGSM_MM[corGSM_MM\$TFRegulator=="YES",1])

      p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, label =  MotifMatrix_matchName, color = TFRegulator)) +
        geom_point() +
        theme_ArchR() +
        geom_vline(xintercept = 0, lty = "dashed") +
        scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
        xlab("Correlation To Gene Score") +
        ylab("Max TF Motif Delta") +
        scale_y_continuous(
          expand = c(0,0),
          limits = c(0, max(corGSM_MM\$maxDelta)*1.05)
        ) +
        geom_text_repel(data = subset(data.frame(corGSM_MM), TFRegulator == "YES"),
                        size = 2,
                        show.legend = FALSE, # otherwise there will be a hidden letter a
                        point.padding = 0.5,
                        force = 200,
                        max.time = 30,
                        segment.size = 0.3,
                        direction = "both",
                        segment.color = "grey50",
                        segment.curvature = -0.1)

      plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)
      },
      error = function(cond) {
        message("TryCatch Error here")
        message(cond)
      }
    )

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
    mkdir -p ./report_jpeg/archr_get_positive_tf_regulator_clusters2
    cp -r ./Plots/jpeg report_jpeg/archr_get_positive_tf_regulator_clusters2 || :
    mkdir ./report_jpeg/archr_get_positive_tf_regulator_clusters2/pdf
    cp ./Plots/*.pdf report_jpeg/archr_get_positive_tf_regulator_clusters2/pdf/ || :


    """
}
