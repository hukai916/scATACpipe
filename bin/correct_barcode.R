#!/usr/bin/env Rscript
# Correct barcode sequence given barcode whitelist and barcode fastq
library(optparse)
option_list = list(
	make_option(c("--barcode_file"), type="character", default=NULL,
							help="full path to barcode fastq file", metavar="character"),
	make_option(c("--whitelist_file"), type="character", default=NULL,
							help="full path to whitelist barcode text file", metavar="character"),
	make_option(c("--reads_per_chunk"), type="integer", default=1000000,
							help="the number of reads to process in one chunk [default=%default]", metavar="character"),
	make_option(c("--path_output_fq"), type="character", default="./",
							help="Output directory [default=%default]", metavar="character"));

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$barcode_file) | is.null(opt$whitelist_file)) {
	print_help(opt_parser)
	stop("At least the barcode_file, whitelist_file must be provided!")
}

library(collections)
library(ShortRead)
library(stringr)
library(purrr)
library(dplyr)
library(campfin)
# library(stringdist) # for two pairs of string vectors, compute a matrix of distance for them, can be fast

# barcode   <- "/Users/kaihu/Projects/workflow/test_data/10x_genomics_10k/atac_hgmm_10k_nextgem_fastqs/atac_hgmm_10k_nextgem_S1_L001_R2_001.fastq.gz"
# whitelist <- "/Users/kaihu/Projects/workflow/test_data/barcodes/737K-cratac-v1.txt.gz"
# reads_per_chunk <- 1000000
# path_output_fq  <- "./"

correct_barcode <- function(barcode_file, whitelist_file, reads_per_chunk, path_output_fq) {
	# Remove preexisting files:
	tem_barcode_corrected_file <- paste0(path_output_fq, "/tem_barcode_corrected_", basename(barcode_file))
	barcode_corrected_file <- paste0(path_output_fq, "/barcode_corrected_", basename(barcode_file))
	file_to_remove <- c(tem_barcode_corrected_file, barcode_corrected_file)
	file.remove(file_to_remove[file.exists(file_to_remove)])


	tem <- readLines(whitelist_file)
	dict_whitelist <- dict(items = as.list(rep(0, length(tem))), keys = as.list(tem))
	dict_invalid	 <- dict()

	# Prepare a dict_1mismatch from whitelist
	## first ensure all whitelist barcode are of the same length (todo)

	total_barcode_records <- countLines(barcode_file) / 4 # faster than countFastq()

	chunk_number <- ceiling(total_barcode_records / reads_per_chunk)
	message(paste0("Total reads: ", total_barcode_records, ", to reduce memory usage, will split into ", chunk_number, " chunk(s), each with a max of ", reads_per_chunk, " reads."))


	invalid_barcode_count <- 0
	f <- FastqStreamer(barcode_file, reads_per_chunk)
	read_chunk <- 1
	while (length(fq <- yield(f))) {
		reads <- fq@sread %>% as.character()
		message(paste0("Processing chunk ", read_chunk, "/", chunk_number, " ..."))
		walk(reads, function(x) {
			if (dict_whitelist$has(x)) {
				dict_whitelist$set(x, dict_whitelist$get(x) + 1)
			} else {
				invalid_barcode_count <<- invalid_barcode_count + 1
				if (!(dict_invalid$has(x))) {
					dict_invalid$set(x, 1)
				} else {
					dict_invalid$set(x, dict_invalid$get(x) + 1)
				}
			}
		})
		invalid_ratio <- format(round(invalid_barcode_count / (read_chunk * reads_per_chunk), 2), nsmall = 2)
		if (invalid_ratio > 0.25) {
			stop("Invalid barcodes exceed 25%, make sure to use the correct barcode whitelist!")
		}
		message(paste0("Processed, invalid barcode: ", invalid_barcode_count, "(", length(dict_invalid$keys()), " unique)/", read_chunk * reads_per_chunk, " (", invalid_ratio, ")"))
		read_chunk <- read_chunk + 1
		# for test only:
		# if (read_chunk == 2) { break }
	}
	close(f)

	keys	 <- unlist(dict_whitelist$keys())
	values <- unlist(dict_whitelist$values())
	tb_whitelist <- tibble(keys, values)

	keys_invalid <- unlist(dict_invalid$keys())
	values			 <- unlist(dict_invalid$values())
	tb_invalid   <- tibble(keys_invalid, values)

	message("For each invalid barcode, find all whitelist barcodes that is 1 mismatch away ...")
	k <- as.character(dict_whitelist$keys())

	wb_wid <- width(k[[1]]) # whitebarcode length
	wb_num <- length(k) # total nubmer of whitebarcodes

	dict_invalid_1mismatch <- dict()

	for (i in seq(wb_wid)) {
		message(paste0("Generating mismatched barcodes at pos ", i, " ..."))
		keys1mismatch 	<- vector(mode = "list", length = 5 * wb_num)
		values1mismatch <- vector(mode = "list", length = 5 * wb_num)
		count <- 0
		for (letter in c("A", "T", "C", "G", "N")) {
			res <- paste0(str_sub(k, 1, i - 1), letter, str_sub(k, i + 1, -1L))
			count <- count + 1
			start <- (count - 1) * wb_num + 1
			end <- count * wb_num
			keys1mismatch[start : end]		<- res
			values1mismatch[start : end]	<- k
		}

		dict_1mismatch <- dict(keys = keys1mismatch, items = values1mismatch)
		message(paste0("  Start processing invalid barcodes ... "))
		walk(keys_invalid, function(x) {
			if (dict_1mismatch$has(x)) {
				if (dict_invalid_1mismatch$has(x)) {
					dict_invalid_1mismatch$set(x, c(dict_1mismatch$get(x), dict_invalid_1mismatch$get(x)))
				} else {
					dict_invalid_1mismatch$set(x, dict_1mismatch$get(x))
				}
			}
		})
		message("Done!")
	}

	message(paste0("Result: ", length(dict_invalid_1mismatch$keys()), "/",  length(keys_invalid), " are 1 mismatch away from whitelist barcodes."))
	message("Determining which whitelist barcode to correct to for each 1 mismatched barcode ...")

	walk(dict_invalid_1mismatch$keys(), function(x) {
		candidates <- dict_invalid_1mismatch$get(x)
		if (!(length(candidates) == 1)) {
			res <- which.max(candidates %>% map(dict_whitelist$get))
			dict_invalid_1mismatch$set(x, candidates[[res]])
		}
	})

	message("Output final barcode corrected fastq ...")

	message(paste0("Total reads: ", total_barcode_records, ", to reduce memory usage, will split into ", chunk_number, " chunk(s), each with a max of ", reads_per_chunk, " reads."))

	f <- FastqStreamer(barcode_file, reads_per_chunk)
	read_chunk <- 1
	total_valid <- 0
	total_1mismatch <- 0
	total_discard <- 0

	while (length(fq <- yield(f))) {
		reads <- fq@sread %>% as.character()
		message(paste0("Processing chunk ", read_chunk, "/", chunk_number, " ..."))
		keep <- vector(mode = "integer", length = reads_per_chunk)

		for (i in seq(reads)) {
			if (dict_whitelist$has(reads[i])) {
				keep[i] <- 1 # indicate barcode perfectly matching whitelist barcode.
			} else if (dict_invalid_1mismatch$has(reads[i])) {
				keep[i]  <- 2 # indicate barcode 1 mismatch away from whitelist barcode.
				reads[i] <- dict_invalid_1mismatch$get(reads[i])
			} else {
				keep[i] <- -1 # indicate barcode more than 2 mismatches from whitelist, should be discarded.
			}
		}

		fq@quality <- fq@quality[which(keep > 0)]
		fq@sread   <- reads[which(keep > 0)] %>% DNAStringSet()
		fq@id      <- fq@id[which(keep > 0)]

		read_chunk <- read_chunk + 1
		valid_count <- sum(keep == 1)
		mismatch1_count <- sum(keep == 2)
		others_count <- sum(keep == -1)

		total_valid <- total_valid + valid_count
		total_1mismatch <- total_1mismatch + mismatch1_count
		total_discard <- total_discard + others_count

		message(paste0("Valid: ", valid_count, "; 1 mismatched: ", mismatch1_count, "; others(discarded): ", others_count))
		writeFastq(fq, paste0(path_output_fq, "/tem_barcode_corrected_", basename(barcode_file)), mode = "a")
	}
	close(f)
	file.rename(paste0(path_output_fq, "/tem_barcode_corrected_", basename(barcode_file)), paste0(path_output_fq, "/barcode_corrected_", basename(barcode_file)))

	message("Barcode correction finished!")
	summary_info <- paste0("Summary (R_correct_barcode): total valid: ", total_valid, "; total corrected (1 mismatch): ", total_1mismatch, "; total discarded: ", total_discard, ".")
	message(summary_info)

	file_connect<-file(paste0("summary_", basename(barcode_file), ".txt"))
	writeLines(summary_info, file_connect)
	close(file_connect)
}

correct_barcode(barcode_file = opt$barcode_file, whitelist_file = opt$whitelist_file, reads_per_chunk = opt$reads_per_chunk, path_output_fq = opt$path_output_fq)
