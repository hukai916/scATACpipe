// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADD_BARCODE_TO_READS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'add_barcode_to_reads', publish_id:'') }
    container "hukai916/sinto_xenial:0.1"

    input:
    // tuple val(sample_name), path(read1_fastq), path(read2_fastq), path(barcode_fastq)
    val sample_name
    path read1_fastq
    path read2_fastq
    path barcode_fastq

    output:
    tuple val(sample_name), path("R1/*barcoded*"), path("R2/*barcoded*"), emit: reads_0

    script:

    """
    # use the first read length from fastq file to determine the length since -b is required by sinto.
    filename=\$(basename -- "$barcode_fastq")
    extension="\${filename##*.}"

    if [[ "\$extension" == "gz" ]]
    then
      barcode_length=\$((zcat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    else
      barcode_length=\$((cat $barcode_fastq || true) | awk 'NR==2 {print length(\$0); exit}')
    fi

    mkdir R1
    cp -P $barcode_fastq R1/
    cp -P $read1_fastq R1/
    sinto barcode $options.args --barcode_fastq R1/$barcode_fastq --read1 R1/$read1_fastq -b \$barcode_length
    rm R1/$barcode_fastq R1/$read1_fastq

    mkdir R2
    cp -P $barcode_fastq R2/
    cp -P $read2_fastq R2/
    sinto barcode $options.args --barcode_fastq R2/$barcode_fastq --read1 R2/$read2_fastq -b \$barcode_length
    rm R2/$barcode_fastq R2/$read2_fastq

    # dev notes:
    # ln \$barcode_fastq R1/ # must be hard link, note hard link won't be created with docker, so use cp instead
    # better to use cp -P to copy soft links
    """
}
