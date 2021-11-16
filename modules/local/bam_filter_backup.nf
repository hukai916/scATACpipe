// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process BAM_FILTER {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'bam_filter', publish_id:'') }

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
    val sample_name
    path bam
    val filter_mitochondrial // if "yes", will filter mitochondrial reads too.

    output:
    val sample_name, emit: sample_name
    path "*.filtered.bam", emit: bam

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes().intdiv(task.cpus)}" : ''
    // Ref: https://gitter.im/nextflow-io/nextflow?at=5a4f8f01ce68c3bc7480d7c5

    // Codes adapted from Haibo Liu.
    if( filter_mitochondrial == 'yes' )
    """
    # Keep only the following reads:
    # 1. Paried reads mapped in the correct orientation.
    # 2. Fragment size ranges from 38 to 2000 bp.
    # 3. The mapq of both reads > 20.
    # 4. Non-mitochondrial reads.

    # Non-mitochondrial chromosome names:
    chromosomes=(`samtools view -H $bam | grep '^@SQ' | perl -n -e 's{.+?SN:([^\t]+).+}{\$1}; if (\$_ ne "MT\n" && \$_ ne "chrM\n") {print}'`)

    # Only output non-mitochondiral reads:
    samtools view -h -b $bam ${chromosomes[@]} | awk 'BEGIN{FS=OFS="\t"} \
    function abs(v) {return v < 0 ? -v : v}; \
    /^@/ || (\$7 == "=" && (\$2 == 99 || \$2 == 147 || \$2 == 83 || \$2 ==163) && abs(\$9) <= 2000 && abs(\$9) >= 38 && \$5 >=20 ) {print}' \
    samtools view -h -b -o ${bam.baseName}.filtered.bam

    """

    else
    """
    # Keep only the following reads:
    # 1. Paried reads mapped in the correct orientation.
    # 2. Fragment size ranges from 38 to 2000 bp.
    # 3. The mapq of both reads > 20.

    samtools view -h -b $bam | awk 'BEGIN{FS=OFS="\t"} \
    function abs(v) {return v < 0 ? -v : v}; \
    /^@/ || (\$7 == "=" && (\$2 == 99 || \$2 == 147 || \$2 == 83 || \$2 ==163) && abs(\$9) <= 2000 && abs(\$9) >= 38 && \$5 >=20 ) {print}' \
    samtools view -h -b -o ${bam.baseName}.filtered.bam

    """
}
