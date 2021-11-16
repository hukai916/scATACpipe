#!/usr/bin/env Rscript

# Create a TxDb .sqlite file given gtf file

library(optparse)
option_list = list(
	make_option(c("--gtf"), type="character", default=NULL,
							help="full path to gtf file", metavar="character"),
	make_option(c("--bsgenome"), type="character", default=NULL,
							help="full path to bsgenome .rds file", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$gtf) | is.null(opt$bsgenome)) {
	print_help(opt_parser)
	stop("The gtf (.gtf.gz or .gtf) and BSgenome (.rds) file must be provided!")
}

library(GenomicFeatures)

# make TxDb from a GTF file:
forgeTxDb <- function(BSgenome, gtf, out_TxDb_dir)
{
    if (missing(BSgenome) | missing(gtf) | missing(out_TxDb_dir))
    {
        stop("All arguments are required: BSgenome, gtf, and out_TxDb_dir!")
    }
    if (!is(BSgenome, "BSgenome")){
        stop(BSgenome, " is not a BSgenome object!")
    }
    if (!file.exists(gtf)){
        stop(gtf, " doesn't exist!")
    }
    if (grepl(".gtf.gz$", gtf))
    {
        in_gtf <- gzfile(gtf, open ="rt")
    } else if (grepl(".gtf$", gtf)) {
        in_gtf <- file(gtf, open = "r")
    } else {
        stop("It seems the GTF file is not a GTF file ",
             "which should with an extension .gtf, or .gtf.gz")
    }
    if (!dir.exists(out_TxDb_dir)){
        dir.create(out_TxDb_dir)
    }

    chrom_len <- seqlengths(BSgenome)
    is_circular <- names(chrom_len) %in% c("chrM", "chrMT", "MT",
                                           "chrPltd", "Pltd")
    chrominfo <- data.frame(chrom = names(chrom_len),
                            length = unname(chrom_len),
                            is_circular = is_circular)
    genome_metadata <- metadata(BSgenome)
    TxDb <- makeTxDbFromGFF(file = in_gtf,
                            format = "gtf",
                            dataSource = NA,
                            organism = NA, # if genome_metadata$genome, unknown organism error
                            taxonomyId = NA,
                            chrominfo = chrominfo,
                            miRBaseBuild = NA)
    close(in_gtf)
    TxDb_file <- file.path(out_TxDb_dir,
                           paste0(genome_metadata$genome, ".TxDb.sqlite"))
    saveDb(TxDb, file = TxDb_file)
    TxDb_file
}

bsgenome <- readRDS(opt$bsgenome)
gtf <- opt$gtf

forgeTxDb(BSgenome = bsgenome, gtf = gtf, out_TxDb_dir = "./")
