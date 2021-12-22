// // Dev notes from preprocess_default:
// include { GET_SAMPLE_NAME_PATH  } from '../modules/local/get_sample_name_path'
// include { GET_SAMPLE_NAME_VAL   } from '../modules/local/get_sample_name_val'
// read1_chunk   = SPLIT_FASTQ.out.read1_fastq.collect().toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()
// read2_chunk   = SPLIT_FASTQ.out.read2_fastq.collect().toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()
// barcode_chunk = SPLIT_FASTQ.out.barcode_fastq.collect().toSortedList( { a, b -> a.name <=> b.name } ).flatten()
// getName() only works for file object, collect()/toSortedList replaces the original filename with the complete order: https://github.com/nextflow-io/nextflow/issues/377
// Here. collect() is a must, otherwise, read1 will be empty when passed to GET_SAMPLE_NAME_PATH: need more reading
// GET_SAMPLE_NAME_PATH (read1_chunk)
// GET_SAMPLE_NAME_VAL (GET_SAMPLE_NAME_PATH.out.sample_name_path)
// sample_name = GET_SAMPLE_NAME_VAL.out.sample_name_val.collect().toSortedList().flatten()
