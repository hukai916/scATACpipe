// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process SPLIT_BAM {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'split_bam', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/sinto_xenial:0.1"
    // cache false

    input:
    // tuple val(sample_name), path(fragment)
    val sample_name
    path tsv // if using path "*.tsv", then in the working dir, the tsv files will be named with numbers, reason unclear.
    path bam
    val barcode_regex

    output:
    path "split_*"
    path "bigWig_*"

    script:

    """
    #bam=rm_dup_${sample_name}.*.bam
    bam=*${sample_name}*.bam
    tsv=(\$(ls *${sample_name}*.tsv))

    for (( i=0; i<\${#tsv[@]}; i++ )); do
      f="\$(basename \${tsv[\$i]} .tsv)"

      # split bam:
      mkdir split_\$f
      cd split_\$f
      samtools index ../\$bam

      if [ "$barcode_regex" == "NA" ]; then
        sinto filterbarcodes -b ../\$bam -c ../\${tsv[\$i]} -p $task.cpus
      else
        sinto filterbarcodes --barcode_regex "$barcode_regex" -b ../\$bam -c ../\${tsv[\$i]} -p $task.cpus
      fi

      # make bw file:
      cd ../
      mkdir bigWig_${sample_name}
      bamfile=(\$(ls -1 split_\$f/*.bam))
      for (( j=0; j<\${#bamfile[@]}; j++ )); do
        filename="\$(basename \${bamfile[\$j]} .bam)"
        samtools index \${bamfile[\$j]}
        bamCoverage -b \${bamfile[\$j]} -o bigWig_${sample_name}/\${filename}.bw -of bigwig $options.bam_coverage 2>bamCoverage.stderr
      done
    done

    """
}
