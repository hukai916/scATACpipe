// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_ARCHRPROJECT_QC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_archrproject_qc', publish_id:'') }
    container "hukai916/r_env_v4.4.1_amd64:0.1"

    input:
    path archr_project
    val archr_thread

    output:
    path archr_project, emit: archr_project
    path "report_jpeg/archr_archrproject_qc", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

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
    mkdir -p ./report_jpeg/archr_archrproject_qc
    cp -r ./Plots/jpeg report_jpeg/archr_archrproject_qc || :
    mkdir ./report_jpeg/archr_archrproject_qc/pdf
    cp ./Plots/*.pdf report_jpeg/archr_archrproject_qc/pdf/ || :

    """
}
