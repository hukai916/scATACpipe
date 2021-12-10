#!/usr/bin/env Rscript
### Use this code to determine the inflection point from 10XGENOMICS scATACseq data.
### This method first truncates out "bumpy" data points from both ends by leveraging the distinct "data point density" pattern. Then, finds the local minimal value of lag diff of read counts.

library(optparse)
option_list = list(
	make_option(c("--freq"), type="character", default=NULL,
							help="Path to file containing read count frequency per cell, must be tsv.", metavar="character"),
	make_option(c("--outfile"), type="character", default="valid_barcode_frequency.txt",
							help="Output filename", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$freq)) {
	print_help(opt_parser)
	stop("The --freq must be provided!")
}
if (file.exists(opt$outfile)) {
  stop("Output file already existS!")
}

# Read in read counts per barcode data and sort by descending order:
freq <- read.csv(opt$freq, sep = "\t", header = FALSE)
freq <- freq[order(-freq$V2),]
x <- log10(c(1:length(freq[, 2])))
y <- log10(freq[, 2])

# Calculate lag diff and get rid of the massive 0s, this won't affect the shape of the curve:
subset0 <- which(diff(y) < 0)

# Visualize the count density to see if there is any pattern that can be leveraged to to remove boundary noises:
# plot(density(x[subset0]))
d <- density(x[subset0])

# Get the truncated range:
subset1 <- which(d$y > max(d$y) - 3 * sd(d$y))
start <- d$x[min(subset1)]
end   <- d$x[max(subset1)]

# Find out valid barcode number and corresponding read counts:
subset2 <- which(x > start & x < end)
num_valid_barcodes <- which(diff(y) == min(diff(y)[subset2])) # number of valid barcodes
cutoff_count <- 10 ** y[diff(y) == min(diff(y)[subset2])] # number of read count that can be used as valid cell cutoff

# Output valid barcodes and their read count frequencies:
write.table(freq[freq$V2 >= cutoff_count,], file = opt$outfile, sep = "\t", quote = F, row.names = F, col.names = F)
