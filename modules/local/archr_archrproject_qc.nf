// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_ARCHRPROJECT_QC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_archrproject_qc', publish_id:'') }

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
    // val archr_filter_ratio
    // val archr_genome
    // val archr_thread

    output:
    path "Plots/TSS-vs-Frags.pdf", emit: pdf_tss_vs_frags
    path "Plots/QC-Sample-Statistics.pdf", emit: pdf_qc_sample_statistics
    path "Plots/QC-Sample-FragSizes-TSSProfile.pdf", emit: pdf_qc_sample_fragsizes_tssprofile
    path "Plots/jpeg", emit: jpeg // to also publish the jpeg folder
    // path "proj_doublet_filtered.rds", emit: archr_project
    path archr_project, emit: archr_project
    // path quality_control, emit: quality_control // if using this syntax, the -resume won't work
    // path "QualityControl", emit: quality_control // using this, the -resume won't work either.
    // This is because the quality_control folder content gets updated after each run, and it will be used as input for itself, so each time, it rerun, the timestamp of this folder is newer.
    path "report_jpeg/archr_archrproject_qc", emit: report
    // path "summary_archrproject_qc.txt", emit: summary

    script:
    // for unknown reason, #!/usr/bin/R + direct R codes won't work
    """
    echo '
    library(ArchR)
    proj <- readRDS("$archr_project", refhook = NULL)

    # Create QC plot: log10(Unique Fragments) vs TSS enrichment score:
    df <- getCellColData(proj, select = c("log10(nFrags)","TSSEnrichment"))
    p <- ggPoint(x = df[,1], y = df[,2], colorDensity = TRUE,
      continuousSet = "sambaNight",
      xlabel = "Log10 Unique Fragments",
      ylabel = "TSS Enrichment",
      xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
      ylim = c(0, quantile(df[,2], probs = 0.99))
      ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
    plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = NULL, addDOC = FALSE)

    # Create QC plot: some statistics: ridge plot, violin plot for TSS enrichment score, violin plot for log10(unique nuclear fragments)
    p1 <- plotGroups(
      ArchRProj = proj,
      groupBy = "Sample",
      colorBy = "cellColData",
      name = "TSSEnrichment",
      plotAs = "ridges"
    )
    p2 <- plotGroups(
      ArchRProj = proj,
      groupBy = "Sample",
      colorBy = "cellColData",
      name = "TSSEnrichment",
      plotAs = "violin",
      alpha = 0.4,
      addBoxPlot = TRUE
    )
    p3 <- plotGroups(
      ArchRProj = proj,
      groupBy = "Sample",
      colorBy = "cellColData",
      name = "log10(nFrags)",
      plotAs = "ridges"
    )
    p4 <- plotGroups(
      ArchRProj = proj,
      groupBy = "Sample",
      colorBy = "cellColData",
      name = "log10(nFrags)",
      plotAs = "violin",
      alpha = 0.4,
      addBoxPlot = TRUE
      )
    plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = NULL, addDOC = FALSE, width = 4, height = 4)

    # Create QC plot: Sample Fragment Size Distribution and TSS Enrichment Profiles
    p5 <- plotFragmentSizes(ArchRProj = proj)
    p6 <- plotTSSEnrichment(ArchRProj = proj)
    plotPDF(p5,p6, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = NULL, addDOC = FALSE, width = 5, height = 5)

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
    mkdir -p report_jpeg/archr_archrproject_qc
    cp -r Plots/jpeg report_jpeg/archr_archrproject_qc
    # cp .command.log summary_archrproject_qc.txt

    """
}
