// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_ARCHRPROJECT_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_archrproject_annotation', publish_id:'') }

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
    path arrowfiles // this will prepare all required files from the list into the working dir
    path gene_annotation
    path genome_annotation

    output:
    // val sample_name, emit: sample_name
    // // path sample_name, emit: output_dir // using this syntax, the -resume won't work
    path "ArchRProject", emit: archrproject_dir
    path "proj.rds", emit: archr_project

    script:
    // for unknown reason, #!/usr/bin/R + direct R codes won't work
    // note the double quote is intentionally used here for echo command, otherwise if singl quote, the $arrow is problematic.
    """
    echo $arrowfiles > arrowfiles.txt
    arrows=\$(cat arrowfiles.txt | sed -e 's/.arrow\s/.arrow", "/g')

    echo "
    library(ArchR)

    genomeAnnotation <- readRDS(\\"$genome_annotation\\")
    geneAnnotation <- readRDS(\\"$gene_annotation\\")

    proj <- ArchRProject(
              ArrowFiles = c(\\"\$arrows\\"),
              geneAnnotation = geneAnnotation,
              genomeAnnotation = genomeAnnotation,
              outputDirectory = \\"ArchRProject\\",
              $options.args
            )

    saveRDS(proj, file = \\"proj.rds\\")
    " > run.R

    Rscript run.R

    """
}
