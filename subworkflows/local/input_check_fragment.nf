//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK_FRAGMENT } from '../../modules/local/samplesheet_check_fragment' addParams( options: params.options )

workflow INPUT_CHECK_FRAGMENT {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK_FRAGMENT ( samplesheet )
        .splitCsv(header: true, sep: ",", strip: true)
        .map { create_fragment_channels(it) }
        .unique()
        .set { fragment }

    emit:
    fragment // channel: [ val(sample_name), path(file_path)]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fragment_channels(LinkedHashMap row) {
  if (!file(row.file_path).exists()) {
      exit 1, "ERROR: Please check input samplesheet -> Fragment file does not exist!\n${row.file_path}"
  }
  array = [ row.sample_name, row.file_path ]

  return array
}
