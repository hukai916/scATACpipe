// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCHR_ARCHRPROJECT_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_archrproject_annotation', publish_id:'') }
    container "hukai916/scatacpipe_downstream:0.2.1"

    input:
    path arrowfiles
    path gene_annotation
    path genome_annotation
    path user_rlib
    val archr_thread

    output:
    path "ArchRProject", emit: archrproject_dir
    path "proj.rds", emit: archr_project

    script:
    // note the double quote is intentionally used here for echo command, otherwise if single quote, the $arrow is problematic.

    """
    echo $arrowfiles > arrowfiles.txt
    arrows=\$(cat arrowfiles.txt | sed -e 's/.arrow\s/.arrow", "/g')

    echo "
    library(ArchR)

    addArchRThreads(threads = $archr_thread)

    # Include the installed custom BSgenome if supplied:
    if (!('$user_rlib' == 'file_token.txt')) {
      .libPaths('user_rlib')
    }

    genomeAnnotation <- readRDS(\\"$genome_annotation\\")
    geneAnnotation <- readRDS(\\"$gene_annotation\\")

    proj <- ArchRProject(
              ArrowFiles = c(\\"\$arrows\\"),
              geneAnnotation = geneAnnotation,
              genomeAnnotation = genomeAnnotation,
              outputDirectory = \\"ArchRProject\\",
              $options.args
            )

    # Fix the gene_symbol containing ENSEMBLE id issue: update codes in create_geneAnnotation_genomeAnnotation.R

    saveRDS(proj, file = \\"proj.rds\\")
    " > run.R

    Rscript run.R

    """
}
