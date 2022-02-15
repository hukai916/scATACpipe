// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_ARCHRPROJECT {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_archrproject', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2"

    input:
    path arrowfiles
    val archr_genome
    val archr_thread

    output:
    path "ArchRProject", emit: archrproject_dir
    path "proj.rds", emit: archr_project
    path "arrowfiles.txt", emit: test_file
    path arrowfiles, emit: arrow_files

    script:
    // note the double quote is intentionally used here for echo command, otherwise if single quote, the $arrow is problematic.

    """
    echo $arrowfiles > arrowfiles.txt
    arrows=\$(cat arrowfiles.txt | sed -e 's/.arrow\s/.arrow", "/g')

    echo "
    library(ArchR)

    addArchRGenome(\\"$archr_genome\\")
    addArchRThreads(threads = $archr_thread)

    proj <- ArchRProject(
              ArrowFiles = c(\\"\$arrows\\"),
              outputDirectory = \\"ArchRProject\\",
              $options.args
            )

    saveRDS(proj, file = \\"proj.rds\\")

    " > run.R

    Rscript run.R

    """
}
