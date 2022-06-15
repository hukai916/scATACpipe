// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CELLRANGER_INDEX {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'cellranger_index', publish_id:'') }
    container "hukai916/cellranger_atat_2.1.0:0.1"

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
