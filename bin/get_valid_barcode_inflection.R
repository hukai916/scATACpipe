#!/usr/bin/env Rscript
### Use this code to determine the inflection point from 10XGENOMICS scATACseq data.
### This method first truncates out "bumpy" data points from both ends by leveraging the distinct "data point density" pattern. Then, finds the local minimal value of lag diff of read counts.

library(optparse)
library(grid)
library(dplyr)

option_list = list(
	make_option(c("--freq"), type="character", default=NULL,
							help="Path to file containing read count frequency per cell, must be tsv.", metavar="character"),
	make_option(c("--outfile"), type="character", default="valid_barcode_frequency.txt",
							help="Output filename", metavar="character"),
	make_option(c("--outplot"), type="character", default="valid_barcode_plot",
							help="Output plotname", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$freq)) {
	print_help(opt_parser)
	stop("The --freq must be provided!")
}
if (file.exists(opt$outfile)) {
  stop("Output file already exists!")
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

# Output a plot showing how valid cell barcodes are determined:
data=data.frame(list(x,y))
names(data) <- c("log10(cell count)", "log10(fragment count)")

min_frag_count <- log10(cutoff_count) # horizontal line
valid_cell_count <- log10(num_valid_barcodes) # vertical line

sp <- ggplot(data=data, aes(x=`log10(cell count)`, y=`log10(barcode frequency per cell)`)) +
				theme(text = element_text(size = 16)) +
				geom_point() +
				geom_vline(xintercept = valid_cell_count, linetype = "dotted", color = "blue", size = 1) +
				geom_hline(yintercept = min_frag_count, linetype = "dotted", color = "red", size = 1)
# add text annotation
annotation_min_frag_count <- paste0("Min fragments for valid cell: ", as.integer(10^min_frag_count), " (10^", format(round(min_frag_count, 2), nsmall = 2), ")")
annotation_valid_cell_count <- paste0("Valid cell count: ", as.integer(10^valid_cell_count), " (10^", format(round(valid_cell_count, 2), nsmall = 2), ")")
grob1 <- grobTree(textGrob(annotation_min_frag_count, hjust = 0, x=0.05, y=(min_frag_count/y[1]), gp=gpar(col="red", fontsize=13, fontface="italic")))
grob2 <- grobTree(textGrob(annotation_valid_cell_count, x=(valid_cell_count/(max(x))), y=0.1, hjust="centre", gp=gpar(col="blue", fontsize=13, fontface="italic")))
sp <- sp + annotation_custom(grob1) + annotation_custom(grob2)
ggsave(paste0(opt$outplot, ".jpeg"),  width = 6, height = 6)
ggsave(paste0(opt$outplot, ".pdf"),  width = 6, height = 6)
