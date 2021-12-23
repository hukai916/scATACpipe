// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTER_BAM {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'filter_bam', publish_id:'') }
    container "hukai916/pysam_xenial:0.1"

    input:
    val sample_name
    path bam
    val filter // if "yes", will filter mitochondrial reads too.

    output:
    tuple val(sample_name), val(chunk_name), path("*.filtered.bam"), emit: sample_name_chunk_name_bam
    val sample_name, emit: sample_name
    path "*.filtered.bam", emit: bam

    script:
    // In case split_fastq is called, chunk_name remains the same as sample_name if not.
    chunk_name = bam.name.split("\\.")[0..-3].join(".") // to document "chunk" info if any

    // Some codes adapted from Haibo Liu, kudos to him!
    if (filter == 'both') // filter out both mitochondiral and "improper reads"
    """
    # Keep only the following reads:
    # 1. Paried reads mapped in the correct orientation.
    # 2. Fragment size ranges from 38 to 2000 bp. (Note, the samtools solution may end up with unpaired reads since we applied other filtering criteria. Therefore, need extract_pair.py.)
    # 3. The mapq of both reads > 20.
    # 4. Non-mitochondrial reads.

    # Non-mitochondrial chromosome names:
    chromosomes=(\$(samtools view -H $bam | grep '^@SQ' | perl -n -e 's{.+?SN:([^\\t]+).+}{\$1}; if (\$_ ne "MT\\n" && \$_ ne "chrM\\n") {print}'))

    # Only output non-mitochondiral reads:
    samtools index $bam
    samtools view -h -@ $task.cpus $options.args $bam \${chromosomes[@]} | awk 'BEGIN{FS=OFS="\\t"} \
    function abs(v) {return v < 0 ? -v : v}; \
    /^@/ || (\$7 == "=" && (\$2 == 99 || \$2 == 147 || \$2 == 83 || \$2 == 163) && abs(\$9) <= 2000 && abs(\$9) >= 38 && \$5 >= 20 ) {print}' \
    > ${bam.baseName}.filtered.tmp.sam

    # Extract paired reads:
    samtools sort -n -@ $task.cpus ${bam.baseName}.filtered.tmp.sam -o ${bam.baseName}.namesrt.tmp.sam
    rm ${bam.baseName}.filtered.tmp.sam
    awk '
        BEGIN { FS=OFS="\\t" }
        FNR == 1 { getline nextline < FILENAME; }
        {
          getline nextline < FILENAME;
          # currentline is \$0, nextline is nextline
          if (\$1 ~ /^@/) { print; next; }
          split(nextline, a);
          if (\$1 == a[1] && \$0 != nextline) {
            print \$0"\\n"nextline;
          }
        }' ${bam.baseName}.namesrt.tmp.sam > ${bam.baseName}.paired.tmp.sam
    rm ${bam.baseName}.namesrt.tmp.sam

    # Output position sorted bam:
    samtools sort -@ $task.cpus ${bam.baseName}.paired.tmp.sam -o ${bam.baseName}.filtered.bam
    rm ${bam.baseName}.paired.tmp.sam

    num_kept=\$(samtools view -c ${bam.baseName}.filtered.bam)
    num_all=\$(samtools view -c $bam)
    num_filtered=\$((num_all - num_kept))

    echo "Summary (bam_filter): total valid: \$num_kept; total filtered: \${num_filtered}." > summary_${bam.baseName}.txt

    """

    else if (filter == 'improper') // filter out only "improper reads"
    """
    # Keep only the following reads:
    # 1. Paired reads mapped in the correct orientation.
    # 2. Fragment size ranges from 38 to 2000 bp.
    # 3. The mapq of both reads > 20.

    samtools index $bam
    samtools view -h -@ $task.cpus $options.args $bam | awk 'BEGIN{FS=OFS="\\t"} \
    function abs(v) {return v < 0 ? -v : v}; \
    /^@/ || (\$7 == "=" && (\$2 == 99 || \$2 == 147 || \$2 == 83 || \$2 == 163) && abs(\$9) <= 2000 && abs(\$9) >= 38 && \$5 >= 20 ) {print}' \
    > ${bam.baseName}.filtered.tmp.sam

    # Extract paired reads:
    samtools sort -n -@ $task.cpus ${bam.baseName}.filtered.tmp.sam -o ${bam.baseName}.namesrt.tmp.sam
    rm ${bam.baseName}.filtered.tmp.sam
    awk '
        BEGIN { FS=OFS="\\t" }
        FNR == 1 { getline nextline < FILENAME; }
        {
          getline nextline < FILENAME;
          # currentline is \$0, nextline is nextline
          if (\$1 ~ /^@/) { print; next; }
          split(nextline, a);
          if (\$1 == a[1] && \$0 != nextline) {
            print \$0"\\n"nextline;
          }
        }' ${bam.baseName}.namesrt.tmp.sam > ${bam.baseName}.paired.tmp.sam
    rm ${bam.baseName}.namesrt.tmp.sam

    # Output position sorted bam:
    samtools sort -@ $task.cpus ${bam.baseName}.paired.tmp.sam -o ${bam.baseName}.filtered.bam
    rm ${bam.baseName}.paired.tmp.sam

    num_kept=\$(samtools view -c ${bam.baseName}.filtered.bam)
    num_all=\$(samtools view -c $bam)
    num_filtered=\$((num_all - num_kept))

    echo "Summary (bam_filter): total valid: \$num_kept; total filtered: \${num_filtered}." > summary_${bam.baseName}.txt

    """

    else // don't apply any filtering
    """
    # Don't apply any default filtering:

    samtools index $bam
    samtools view -h -b -@ $task.cpus $options.args $bam -o ${bam.baseName}.filtered.bam

    num_kept=\$(samtools view -c ${bam.baseName}.filtered.bam)
    num_all=\$(samtools view -c $bam)
    num_filtered=\$((num_all - num_kept))

    echo "Summary (bam_filter): total valid: \$num_kept; total filtered: \${num_filtered}." > summary_${bam.baseName}.txt

    """

}
