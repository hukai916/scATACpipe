// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process PREP_GENOME {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'prep_genome', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/bwa_xenial:0.1"
    // cache false

    input:
    path genome_fasta
    val genome_name

    output:
    path "Genome.primary.chrPrefixed.fa.gz", emit: genome_fasta
    val genome_name, emit: genome_name

    """
    if [[ $genome_fasta == *.gz ]]; then
      gunzip -c $genome_fasta > genome.fa
    else
      if [[ $genome_fasta != genome.fa ]]; then
        mv $genome_fasta genome.fa
      fi
    fi

    samtools faidx genome.fa

    # get primary genome:
    cat genome.fa.fai | cut -f 1 | grep -v -P "_alt\$|_fix\$|_hap\\d*\$" | xargs -n 1 -I {} samtools faidx genome.fa {} >> primary_genome.fa

    # add 'chr' to each chromosome:
    cat primary_genome.fa | perl -p -e 's{>(chr)?([^\\s+]+).*}{>chr\$2} if /^>/' | gzip > Genome.primary.chrPrefixed.fa.gz

    rm genome.fa genome.fa.fai primary_genome.fa
    """
}
