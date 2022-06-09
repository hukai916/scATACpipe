// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_PAIRWISE_TEST_CLUSTERS2 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_pairwise_test_clusters2', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

    input:
    path archr_project
    val archr_thread

    output:
    path "markerTest.rds", emit: archr_marker_test
    path "test_group.txt", emit: test_group
    path "report_jpeg/archr_pairwise_test_clusters2", emit: report

    script:

    """
    echo '
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    proj <- readRDS("$archr_project")

    # Determine the useGroups and bgdGroups, default to the first and second:
    if ("$options.use_groups" == "default") {
      useGroups <- unique(proj\$Clusters2)[1]
    } else {
      useGroups <- "$options.use_groups"
    }
    if ("$options.bgd_groups" == "default") {
      bgdGroups <- unique(proj\$Clusters2)[2]
    } else {
      bgdGroups <- "$options.bgd_groups"
    }
    sink(file = "test_group.txt") # first line
    cat(useGroups, "\n")
    cat(bgdGroups, "\n")
    sink()

    markerTest <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = "Clusters2",
      useGroups = useGroups,
      bgdGroups = bgdGroups,
      $options.args
    )
    saveRDS(markerTest, file = "markerTest.rds")

    pma <- markerPlot(seMarker = markerTest, name = useGroups, plotAs = "MA", $options.cutoff)
    pv <- markerPlot(seMarker = markerTest, name = useGroups, plotAs = "Volcano", $options.cutoff)

    plotPDF(pma, pv, name = paste0(useGroups, "-vs-", bgdGroups, "-Markers-MA-Volcano"), width = 5, height = 5, ArchRProj = NULL, addDOC = FALSE)

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
    mkdir -p ./report_jpeg/archr_pairwise_test_clusters2
    cp -r ./Plots/jpeg report_jpeg/archr_pairwise_test_clusters2
    mkdir ./report_jpeg/archr_pairwise_test_clusters2/pdf
    cp ./Plots/*.pdf report_jpeg/archr_pairwise_test_clusters2/pdf/


    """
}
