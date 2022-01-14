// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AMULET_DETECT_DOUBLETS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'amulet_detect_doublets', publish_id:'') }
    container "hukai916/amulet_xenial:0.1"

    input:
    tuple val(sample_name), path(fragment)
    path amulet_rmsk_bed
    path amulet_autosomes

    output:
    path "cells_filter_${sample_name}.txt", emit: cells_filter

    script:

    """
    # prepare singlecell.csv that is required by AMULET:
    # since fragments are already filtered, all fragments are from valid cells.
    echo "barcode,is__cell_barcode" > singlecell_${sample_name}.csv
    zcat $fragment | grep -v "#" | cut -f 4 | sort | uniq | awk '{ print \$1",1" }' >> singlecell_${sample_name}.csv

    mkdir doublets_${sample_name}

    # find overlappings:
    FragmentFileOverlapCounter.py $options.args $fragment singlecell_${sample_name}.csv $amulet_autosomes doublets_${sample_name}
    # detect multiplets :
    AMULET.py --rfilter $amulet_rmsk_bed doublets_${sample_name}/Overlaps.txt doublets_${sample_name}/OverlapSummary.txt doublets_${sample_name}
    # prepare a ArchR compatible cellsFilter object
    cat doublets_${sample_name}/MultipletBarcodes_01.txt | awk -v sample_name="${sample_name}" '{ print sample_name"#"\$0}' > cells_filter_${sample_name}.txt

    """
}
