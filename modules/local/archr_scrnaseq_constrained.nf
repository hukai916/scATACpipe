// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_SCRNASEQ_CONSTRAINED {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_scrnaseq_constrained', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

    input:
    path archr_project
    path obj_scrnaseq
    path archr
    val group_list
    val archr_thread

    output:
    path "proj_scrnaseq_constrained.rds", emit: archr_project
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    path "report_jpeg/archr_scrnaseq_constrained", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(collections)

    addArchRThreads(threads = $archr_thread)

    seRNA <- readRDS("$obj_scrnaseq")
    proj <- readRDS("$archr_project", refhook = NULL)

    # Prepare for the groupList from params.archr_scrnaseq_grouplist:
    cM <- as.matrix(confusionMatrix(proj\$Clusters, proj\$predictedGroup_Un))
    preClust <- colnames(cM)[apply(cM, 1 , which.max)]
    dict1 = dict(list($group_list))

    for (x in dict1\$keys()) {
      c_name <- paste0("c_", x)
      assign(c_name, paste0(dict1\$get(x), collapse = "|"))

      cluster_name <- paste0("cluster_", x)
      assign(cluster_name, rownames(cM)[grep(get(c_name), preClust)])

      rna_name <- paste0("rna_", x)
      assign(rna_name, colnames(seRNA)[grep(get(c_name), colData(seRNA)\$BioClassification)])
    }

    groupList <- SimpleList()

    for (x in dict1\$keys()) {
      c_name	       <- paste0("c_", x)
      cluster_name   <- paste0("cluster_", x)
      rna_name			 <- paste0("rna_", x)
      group_name		 <- paste0("group_", x)

      assign(group_name, SimpleList(
        ATAC = proj\$cellNames[proj\$Clusters %in% get(cluster_name)],
        RNA = get(rna_name)
    	))

      if (length(groupList) == 0) {
        #groupList <- SimpleList(assign(group_name, get(group_name)))
        #groupList <- SimpleList(name1 = get(group_name))
        groupList <- eval(parse(text=paste0("SimpleList(", eval(group_name), " = get(group_name))")))
      } else {
        #tem <- SimpleList(assign(group_name, get(group_name)))
        #tem <- SimpleList(eval(parse(text=paste0(eval(group_name), " = get(group_name)"))))
        tem <- eval(parse(text=paste0("SimpleList(", eval(group_name), " = get(group_name))")))
        # Reference: https://stackoverflow.com/questions/1743698/evaluate-expression-given-as-a-string
        groupList <- append(groupList, tem)
      }
    }

    # Perform constrained integration
    proj2 <- addGeneIntegrationMatrix(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      matrixName = "GeneIntegrationMatrix",
      reducedDims = "IterativeLSI",
      seRNA = seRNA,
      addToArrow = FALSE,
      groupList = groupList,
      groupRNA = "BioClassification",
      nameCell = "predictedCell_Co",
      nameGroup = "predictedGroup_Co",
      nameScore = "predictedScore_Co"
    )
    proj2 <- addGeneIntegrationMatrix(
      ArchRProj = proj2,
      useMatrix = "GeneScoreMatrix",
      matrixName = "GeneIntegrationMatrix",
      reducedDims = "IterativeLSI",
      seRNA = seRNA,
      addToArrow = TRUE,
      force= TRUE,
      groupList = groupList,
      groupRNA = "BioClassification",
      nameCell = "predictedCell",
      nameGroup = "predictedGroup",
      nameScore = "predictedScore"
    )

    # Add impute weights
    proj2 <- addImputeWeights(proj2)

    # markerGenes: default to use first 3 marker_genes inferred from "Clusters"
    if (!("$options.marker_genes" == "default")) {
      markerGenes <- str_trim(str_split("$options.marker_genes", ",")[[1]], side = "both")
    } else {
      markersGS <- getMarkerFeatures(
                     ArchRProj = proj,
                     groupBy = "Clusters",
                     useMatrix = "GeneScoreMatrix",
                     bias = c("TSSEnrichment", "log10(nFrags)"),
                     testMethod = "wilcoxon"
                   )
      markerList <- getMarkers(markersGS)

      markerGenes <- c()
      for (cluster in markerList@listData) {
        markerGenes <- c(markerGenes, cluster\$name)
      }
      sel <- min(length(markerGenes), 3)
      markerGenes <- markerGenes[1:sel]
    }
    # markerGenes must also be a subset of .getFeatureDF(getArrowFiles(proj2), "GeneIntegrationMatrix")\$name
    devtools::load_all("ArchR") # to use .getFeatureDF() function
    geneDF <- .getFeatureDF(getArrowFiles(proj2), "GeneScoreMatrix")
    geneDF <- geneDF[geneDF\$name %in% rownames(seRNA), , drop = FALSE]
    markerGenes <- markerGenes[markerGenes %in% geneDF\$name]
    devtools::unload("ArchR")
    library(ArchR)
    # note that after load_all, the proj will not be correctly recognized as S4 object, have to unload and re-library ArchR.

    if ("Harmony" %in% names(proj@reducedDims)) {
      reducedDims <- "Harmony"
    } else if ("IterativeLSI" %in% names(proj@reducedDims)) {
      reducedDims <- "IterativeLSI"
    }

    embedding <- paste0("UMAP_", reducedDims)

    p1 <- plotEmbedding(
            ArchRProj = proj2,
            embedding = embedding,
            colorBy = "GeneIntegrationMatrix",
            name = markerGenes,
            continuousSet = "horizonExtra",
            imputeWeights = getImputeWeights(proj2)
          )
    plotPDF(plotList = p1,
              name = paste0("Plot-", embedding, "-Marker-Genes-RNA-W-Imputation.pdf"),
              ArchRProj = NULL,
              addDOC = FALSE, width = 5, height = 5
            )

    # Label scATAC-seq clusters with scRNA-seq information
    cM <- confusionMatrix(proj2\$Clusters, proj2\$predictedGroup)
    labelOld <- rownames(cM)
    labelNew <- colnames(cM)[apply(cM, 1, which.max)]
    proj2\$Clusters2 <- mapLabels(proj2\$Clusters, newLabels = labelNew, oldLabels = labelOld)
    p1 <- plotEmbedding(
            ArchRProj = proj2,
            embedding = embedding,
            colorBy = "cellColData",
            name = "Clusters2"
          )
    plotPDF(p1, name = paste0("Plot-", embedding, "-Remap-Clusters.pdf"), ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

    saveRDS(proj2, file = "proj_scrnaseq_constrained.rds")

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
    mkdir -p report_jpeg/archr_scrnaseq_constrained
    cp -r Plots/jpeg report_jpeg/archr_scrnaseq_constrained

    """
}
