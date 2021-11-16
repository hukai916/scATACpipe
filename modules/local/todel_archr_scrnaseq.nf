// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_SCRNASEQ {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_scrnaseq', publish_id:'') }

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
    path obj_scrnaseq

    output:
    path "proj_scrnaseq.rds", emit: archr_project
    path "Plots/Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", emit: pdf_umap_marker_genes_rna_w_imputation
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder

    script:

    """
    echo '
    library(ArchR)
    seRNA <- readRDS("$obj_scrnaseq")
    proj <- readRDS("$archr_project", refhook = NULL)

    # Perform unconstrained integration
    proj2 <- addGeneIntegrationMatrix(
        ArchRProj = proj,
        seRNA = seRNA,
        nameCell = "predictedCell",
        nameGroup = "predictedGroup",
        nameScore = "predictedScore",
        $options.args
    )

    # Add impute weights
    proj2 <- addImputeWeights(proj2)
    saveRDS(proj2, file = "proj_scrnaseq.rds")

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

    """
}
