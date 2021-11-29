// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_SAMPLE_NAME_PATH {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'get_sample_name_path', publish_id:'') }
    container "hukai916/sinto_xenial:0.1"

    input:
    path sample

    output:
    path "*.token", emit: sample_name_path


    script:

    """
    #!/usr/bin/env python

    sample_name = "_".join("$sample".split("_")[1:-2])

    with open(sample_name + ".token", "w") as f:
      pass

    """
}
