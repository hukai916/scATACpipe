// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CELLRANGER_ATAC_COUNT {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'cellranger_count', publish_id:'') }
    container "hukai916/cellranger_atat_2.0.0:0.1"

    input:
    val sample_name
    path sample_files
    path reference

    output:
    val sample_name, emit: sample_name
    path "cellranger_atac_count_*/outs/fragments.tsv.gz", emit: fragments
    tuple val(sample_name), path("cellranger_atac_count_*/outs/fragments.tsv.gz"), emit: ch_fragment
    path "cellranger_atac_count_*", emit: cellranger_atac_count
    path "cellranger_atac_count_*/outs/*_possorted_bam.bam", emit: bam

    script:
    def avail_mem = task.memory ? "${ (task.memory.toBytes().intdiv(1073741824) * 0.9).toInteger() }" : ''

    """

    # the fastq file name must not contain special characters other than dash, underscore, digit; dot is not allowed
    # the --id must not contain dot either:

    # prepare for output folder:
    outfolder=cellranger_atac_count_\$( echo $sample_name | tr '.' '_' ) # just in case

    # prepare for input fastq folder:
    infastq=input_fastq_\$outfolder
    mkdir \$infastq
    cp ${sample_name}_S1_L*_*_001.fastq.gz \$infastq/

    # cellranger_atac is very strict regarding fastq nomenclature, need to match it
    # note that the Lane number may not match with the original, but it should not matter
    cd \$infastq/
    # get number of lanes
    lanes=(\$(ls -d *_R1_001.fastq.gz))
    lanes=\${#lanes[@]}

    # get unique sample_count
    sample_count=array()
    sample_fastq=(\$(ls -d *.fastq.gz))
    for fastq in "\${sample_fastq[@]}"
    do
      [[ \$fastq =~ _L([0-9]+)_R[0-9]_(001.fastq.gz) ]]
      sample_count+=(\${BASH_REMATCH[1]})
    done
    uniq_sample_count=(\$(echo "\${sample_count[@]}" | tr ' ' '\\n' | sort -u | tr '\\n' ' '))

    # rename sample_count to use formatted lane number
    for sample_count in "\${sample_count[@]}"
    do
      printf -v lane "%03d" \$sample_count
      mv ${sample_name}_S1_L\${sample_count}_R1_001.fastq.gz ${sample_name}_S1_L\${lane}_R1_001.fastq.gz
      mv ${sample_name}_S1_L\${sample_count}_R2_001.fastq.gz ${sample_name}_S1_L\${lane}_R2_001.fastq.gz
      mv ${sample_name}_S1_L\${sample_count}_R3_001.fastq.gz ${sample_name}_S1_L\${lane}_R3_001.fastq.gz
    done

    cellranger-atac count $options.args \
    --id \$outfolder \
    --fastqs \$infastq \
    --reference $reference \
    --localcores $task.cpus \
    --localmem $avail_mem

    # rename the output bam file for split_bam module:
    mv \${outfolder}/outs/possorted_bam.bam \${outfolder}/outs/outs/${sample_name}_possorted_bam.bam
    mv \${outfolder}/outs/possorted_bam.bam.bai \${outfolder}/outs/outs/${sample_name}_possorted_bam.bam.bai

    """
}
