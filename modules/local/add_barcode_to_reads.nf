// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ADD_BARCODE_TO_READS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'add_barcode_to_reads', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }
    // container "hukai916/add_barcode_to_reads:0.1"
    container "hukai916/sinto_xenial:0.1"

    // cache false

    input:
    val sample_name
    path barcode_fastq
    path read1_fastq
    path read2_fastq

    output:
    val sample_name, emit: sample_name
    path "R1/*barcoded*", emit: read1_fastq
    path "R2/*barcoded*", emit: read2_fastq

    script:

    """
    # use the first read length from fastq file to determine the length since -b is required by sinto.
    filename=\$(basename -- "$barcode_fastq")
    extension="\${filename##*.}"

    if [[ "\$extension" == "gz" ]]
    then
      # barcode_length=\$(zcat $barcode_fastq | awk '{if(NR==2) print length(\$1)}')
      # Below is more efficient: but also exit 141 when running on NF.
      # barcode_length=\$(zcat $barcode_fastq | awk 'NR==2 { print length(\$1); exit }')
      # Below works too, but not with head -n 1, since head breaks the pipe and exit with 141.
      barcode_length=\$(zcat $barcode_fastq | awk '{if(NR%4==2) print length(\$1)}' | tail -n 1)
    else
      barcode_length=\$(cat $barcode_fastq | awk 'NR==2 { print length(\$1); exit }')
    fi

    mkdir R1
    # ln $barcode_fastq R1/ # must be hard link, note hard link won't be created with docker, so use cp instead
    # ln $read1_fastq R1/
    cp $barcode_fastq R1/
    cp $read1_fastq R1/
    sinto barcode $options.args --barcode_fastq R1/$barcode_fastq --read1 R1/$read1_fastq -b \$barcode_length
    rm R1/$barcode_fastq R1/$read1_fastq

    mkdir R2
    # ln $barcode_fastq R2/ # must be hard link
    # ln $read2_fastq R2/
    cp $barcode_fastq R2/
    cp $read2_fastq R2/
    sinto barcode $options.args --barcode_fastq R2/$barcode_fastq --read1 R2/$read2_fastq -b \$barcode_length
    rm R2/$read2_fastq R2/$barcode_fastq
    """
}
