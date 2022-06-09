// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_FOOTPRINTING_CLUSTERS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_footprinting_clusters', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

    input:
    path archr_project
    path user_rlib
    val archr_thread

    output:
    path "report_jpeg/archr_footprinting_clusters", emit: report

    script:

    """
    echo '
    library(ArchR)
    library(stringr)
    .libPaths("user_rlib") # for user installed packages

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project", refhook = NULL)
    dir.create("Plots")
    dir.create("Plots/jpeg")

    # Footprinting of motif:
    motifPositions <- getPositions(proj)

    if ("$options.motifs" == "default") {
      plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
      VarDev     <- getVarDeviations(proj, name = "MotifMatrix", plot = FALSE)
      motifs     <- VarDev\$name[1:min(3, length(VarDev\$name))]
    } else {
      motifs <- str_trim(str_split("$options.motifs", ",")[[1]], side = "both")
    }
    markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

    seFoot <- getFootprints(
                ArchRProj = proj,
                positions = motifPositions[markerMotifs],
                groupBy   = "Clusters"
              )
    plotName <- paste0("Footprints", "-", "$options.norm_method", "-Bias")
    out <- tryCatch(
      expr = {
        plotFootprints(
          seFoot      = seFoot,
          ArchRProj   = proj,
          normMethod  = "$options.norm_method",
          plotName    = plotName,
          $options.args
        )

        file.copy(from = paste0(getOutputDirectory(proj), "./Plots/", "Footprints", "-", "$options.norm_method", "-Bias", ".pdf"), to = paste0("./Plots/", "Footprints", "-", "$options.norm_method", "-Bias", ".pdf"))
      },
      error = function(e) {
        return("Footprint plotting failed.")
      }
    )

    # Footprinting of TSS (custom) Features
    seTSS <- getFootprints(
              ArchRProj = proj,
              positions = GRangesList(TSS = getTSS(proj)),
              groupBy   = "Clusters",
              flank     = $options.tss_flank
             )
    out <- tryCatch(
      expr = {
        plotFootprints(
          seFoot      = seTSS,
          ArchRProj   = proj,
          normMethod  = "$options.tss_norm_method",
          plotName    = paste0("TSS-", "$options.tss_norm_method", "-Normalization"),
          addDOC      = FALSE,
          flank       = $options.tss_flank,
          flankNorm   = $options.flank_norm
          )

          file.copy(from = paste0(getOutputDirectory(proj), "./Plots/", "TSS-", "$options.tss_norm_method", "-Normalization", ".pdf"), to = paste0("./Plots/", "TSS-", "$options.tss_norm_method", "-Normalization", ".pdf"))
      },
      error = function(e) {
        return("Footprint plotting failed.")
      }
    )

    ' > run.R

    Rscript run.R

    # Convert to jpeg:
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
    mkdir -p ./report_jpeg/archr_footprinting_clusters
    cp -r ./Plots/jpeg report_jpeg/archr_footprinting_clusters
    mkdir ./report_jpeg/archr_footprinting_clusters/pdf
    cp ./Plots/*.pdf report_jpeg/archr_footprinting_clusters/pdf

    """
}
