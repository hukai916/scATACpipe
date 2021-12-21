// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MATCH_CHUNK {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'match_chunk', publish_id:'') }
    container "hukai916/sinto_xenial:0.1"

    input:
    path read1_chunk
    path read2_chunks
    path barcode_chunks

    output:
    tuple val("${sample_name}"), path(read1_chunk), path("${read2_chunk}"), path("${barcode_chunk}"), emit: chunk

    script:
    sample_name = read1_chunk.name.split("_")[1..-3].join("_")
    suffix      = read1_chunk.name.split("_")[-2..-1].join("_")
    read2_chunk = "R2_" + sample_name + "_" + suffix
    barcode_chunk = "barcode_" + sample_name + "_" + suffix

    """
    echo $sample_name > result.txt
    echo $read2_chunk >> result.txt

    """

}
