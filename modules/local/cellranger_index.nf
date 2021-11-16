// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process CELLRANGER_INDEX {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'cellranger_index', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/cellranger_atat_2.0.0:0.1"
    // cache false

    input:
    path genome_fasta
    path gtf
    val genome_name

    output:
    path "genome_index", emit: index_folder

    script:

    """
    # Unzip genome_fasta and gtf file
    gunzip -c $genome_fasta > genome.fa
    gunzip -c $gtf > annotation.gtf

    # Prepare config file:
    mem=\$(echo \"$task.memory\" | sed "s/ GB//")
    echo '{' >> index.config
    echo '    organism: \"$genome_name\"' >> index.config
    echo '    genome: [ \"genome_index\"]' >> index.config
    echo '    input_fasta: [\"genome.fa\"]' >> index.config
    echo '    input_gtf: [\"annotation.gtf\"]' >> index.config
    #echo '    non_nuclear_contigs: [\"MT\"]' >> index.config
    echo 'nthreads: $task.cpus' >> index.config
    echo "memgb: \$mem" >> index.config

    echo '}' >> index.config

    # In case GTF contains config names that are not in the genome.fa:
    

    # Make ref:
    cellranger-atac mkref --config=index.config

    """
}
