// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_PEAK2GENELINKAGE_CLUSTERS2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_peak2genelinkage_clusters2', publish_id:'') }
    container "hukai916/r_sc:0.5"

    input:
    path archr_project
    val archr_thread

    output:
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "archr_peak2genelinkage_clusters2", emit: res_dir
    path "report_jpeg/archr_peak2genelinkage_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)

    proj2 <- addPeak2GeneLinks(
      ArchRProj = proj,
      reducedDims = "IterativeLSI"
      )

    p2g <- getPeak2GeneLinks(
        ArchRProj = proj2,
        returnLoops = TRUE,
        $options.args
    )

    # Draw heatmap: default to use first 10 marker_genes
    if (!("$options.marker_genes" == "default")) {
      markerGenes <- str_trim(str_split("$options.marker_genes", ",")[[1]], side = "both")
      markerGenes <- unique(markerGenes)
    } else {
      markerGenes <- c()
      for (cluster in markerList@listData) {
        markerGenes <- c(markerGenes, cluster\$name)
      }
      markerGenes <- unique(markerGenes)
      sel <- min(length(markerGenes), 10)
      markerGenes <- markerGenes[1:sel]
    }

    markerGenes_clean <- markerGenes

    all_id <- getGenes(proj)\$gene_id
    all_symbol <- getGenes(proj)\$symbol
    all_symbol_cleaned <- character(length(all_id))
    for (i in 1:length(all_id)) {
      all_symbol_cleaned[i] <- str_remove(all_symbol[i], paste0("_", all_id[i]))
      markerGenes_clean <- str_remove(markerGenes_clean, paste0("_", all_id[i])) # not very efficient, but works
    }

    markerGenes2labeled <- sort(markerGenes_clean[markerGenes_clean %in% all_symbol_cleaned])
    markerGenes_raw <- sort(all_symbol[all_symbol_cleaned %in% markerGenes_clean])

    all_symbol_cleaned_unique <- unique(all_symbol_cleaned)
    all_symbol_unique <- all_symbol[match(all_symbol_cleaned_unique, all_symbol_cleaned)]
    markerGenes_raw <- sort(all_symbol_unique[all_symbol_cleaned_unique %in% markerGenes_clean])




    markerGenes <- c($options.marker_genes)
    p <- plotBrowserTrack(
      ArchRProj = proj2,
      groupBy = "Clusters2",
      geneSymbol = markerGenes,
      upstream = 50000,
      downstream = 50000,
      loops = getPeak2GeneLinks(proj2)
      )
    plotPDF(plotList = p,
      name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",
      ArchRProj = NULL,
      addDOC = FALSE, width = 5, height = 5)

    p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Clusters2")
    plotPDF(p,
      name = "Heatmap-Marker-Genes-with-Peak2GeneLinks.pdf",
      ArchRProj = NULL,
      addDOC = FALSE, width = 5, height = 12)

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

    # Copy to res_dir:
    mkdir archr_peak2genelinkage_clusters2
    cp -r Plots archr_peak2genelinkage_clusters2/

    # For reporting:
    mkdir -p report_jpeg/archr_peak2genelinkage_clusters2
    cp -r archr_peak2genelinkage_clusters2 report_jpeg/archr_peak2genelinkage_clusters2

    """
}
