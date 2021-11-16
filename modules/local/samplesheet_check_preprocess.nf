// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Reformat design file and check validity
 */
process SAMPLESHEET_CHECK_PREPROCESS {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    container "hukai916/miniconda3_bio:0.1"

    input:
    path samplesheet

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/scatacseqflow/bin/
    """
    check_samplesheet_preprocess.py $samplesheet samplesheet.valid.csv
    """
}
