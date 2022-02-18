// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLIT_BAM {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'split_bam', publish_id:'') }

    container "hukai916/sinto_xenial:0.1"

    input:
    val sample_name
    path tsv // if using path "*.tsv", then in the working dir, the tsv files will be named with numbers, reason unclear.
    path bam

    output:
    path "split_*"
    path "bigWig_*"

    script:

    """
    bam=*${sample_name}*.bam
    tsv=(\$(ls *${sample_name}*.tsv))

    for (( i=0; i<\${#tsv[@]}; i++ )); do
      f="\$(basename \${tsv[\$i]} .tsv)"

      # split bam:
      mkdir split_\$f
      cd split_\$f
      samtools index ../\$bam

      if [[ "$options.barcode_regex" == "" || "$options.barcode_regex" == "null" || "$options.barcode_regex" == "false" ]]; then
        sinto filterbarcodes -b ../\$bam -c ../\${tsv[\$i]} -p $task.cpus
      else
        sinto filterbarcodes --barcode_regex "$options.barcode_regex" -b ../\$bam -c ../\${tsv[\$i]} -p $task.cpus
      fi

      # make bw file:
      cd ../
      mkdir bigWig_\$f
      bamfile=(\$(ls -1 split_\$f/*.bam))
      for (( j=0; j<\${#bamfile[@]}; j++ )); do
        filename="\$(basename \${bamfile[\$j]} .bam)"
        samtools index \${bamfile[\$j]}
        bamCoverage -b \${bamfile[\$j]} -o bigWig_\$f/\${filename}.bw -of bigwig $options.bam_coverage 2>bamCoverage.stderr
      done
    done

    """
}
