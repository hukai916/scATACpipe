//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK_PREPROCESS } from '../../modules/local/samplesheet_check_preprocess' addParams( options: params.options )

workflow INPUT_CHECK_PREPROCESS {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK_PREPROCESS ( samplesheet )
        .splitCsv(header: true, sep: ",", strip: true)
        .map {
          row ->
            [ row.sample_name, row.path_fastq_1, row.path_fastq_2, row.path_barcode ]
        }
        .unique()
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}
