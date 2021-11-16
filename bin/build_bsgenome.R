#!/usr/bin/env Rscript

# Create a BSgenome .rds file given genome fasta and genome name
library(optparse)
option_list = list(
	make_option(c("--genome_fasta"), type="character", default=NULL,
							help="full path to genome fasta", metavar="character"),
	make_option(c("--genome_name"), type="character", default="Custom genome",
							help="Latin name of the genome", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$genome_fasta) | is.null(opt$genome_name)) {
	print_help(opt_parser)
	stop("At least the genome_fasta must be provided!")
}

library(BSgenome)
library(devtools)

generate_multifasta <- function(path_single_fasta_genome, out_dir)
{
    if (missing(path_single_fasta_genome) | missing(out_dir))
    {
        stop("path_single_fasta_genome and out_dir are required!")
    }
    if (!file.exists(path_single_fasta_genome)){
        stop(path_single_fasta_genome, " doesn't exists!")
    }
    if (!dir.exists(out_dir))
    {
        dir.create(out_dir, recursive = TRUE)
    }

    if (grepl(".(fa|fasta).gz$", path_single_fasta_genome))
    {
        in_fasta <- gzfile(path_single_fasta_genome, open ="rt")
    } else if (grepl(".(fa|fasta)$", path_single_fasta_genome)) {
        in_fasta <- file(path_single_fasta_genome, open = "r")
    } else {
        stop("It seems the genome sequence file is not a fasta file ",
             "which should with an extension .fa, .fasta, .fa.gz or .fasta.gz")
    }

    f <- ""
    while (length({line <- readLines(in_fasta, n = 1, warn = FALSE)}) > 0)
    {
        line <- trimws(line)
        if (grepl("^>", line))
        {
            f <- gzfile(file.path(out_dir,
                                  gsub("^>(chr)?([^\\s]+)", "chr\\2.fa.gz",
                                       line, perl = TRUE)),
                        "w")
            writeLines(gsub("^>(chr)?([^\\s]+)",">chr\\2", line, perl = TRUE), f)
        } else {
            writeLines(line, f)
        }
    }
    close(f)
    close(in_fasta)
    out_dir
}

## generate a seed file for building a BSgenome package
generate_seed_file <- function(path_to_multifasta,
                               latin_name,
                               common_name,
                               genome_build,
                               seed_file_name,
                               fasta_url,
                               source = "Ensembl",
                               version = "1.0.0")
{
    if (missing(path_to_multifasta) | missing(latin_name) |
        missing(common_name) | missing(genome_build) |
        missing(seed_file_name) | missing(fasta_url))
    {
        stop("All arguments except source and version are required!")
    }
    if (!dir.exists(path_to_multifasta)) {
        stop("Path to multifasta ", path_to_multifasta, " doesn't exist!")
    }
    chr_fa_files <- dir(path_to_multifasta, ".fa.gz$")
    if (length(chr_fa_files) < 1)
    {
        stop("There is no multiple fasta files in the directory ",
             path_to_multifasta, "!")
    }
    sink(seed_file_name)
    BSgenomeObjname <- gsub("^(.).*\\s+(.+)", "\\1\\2", latin_name)
    package_name <- paste("BSgenome", BSgenomeObjname, source, genome_build, sep = ".")
    cat(paste0("Package: ", package_name, "\n"))
    cat(paste("Title: Full genome sequences for", latin_name,
             paste0("(", source), "version", paste0(genome_build,")\n")))
    cat(paste("Description: Full genome sequences for", latin_name,
             paste0("(", common_name, ")"), "as provided by",
             source, paste0("(", genome_build, ")"),
             "and stored in Biostrings objects.\n"))
    cat(paste0("Version: ", version, "\n"))
    cat(paste0("organism: ", latin_name, "\n"))
    cat(paste0("common_name: ", common_name,"\n"))
    cat(paste0("provider: ", source, "\n"))
    cat("release_date: May, 2007\n")
    cat(paste0("genome: ", genome_build, "\n"))
    cat(paste0("source_url: ", fasta_url, "\n"))
    cat(paste0("BSgenomeObjname: ", BSgenomeObjname, "\n"))
    cat(paste0("organism_biocview: ",
               gsub("\\s+", "_", latin_name, perl = TRUE), "\n"))
    chromosome_names <- gsub(".fa.gz$", "", dir(path_to_multifasta, "fa.gz$"))
    cat(paste0('seqnames: c("', paste(chromosome_names, collapse = '","'), '")\n'))

    circ_seqs <- c("chrM", "MT", "Pltd", "chrPltd")
    circ_seqs <- circ_seqs[circ_seqs %in% chromosome_names]
    cat(paste0('circ_seqs: c("', paste(circ_seqs, collapse = '","'), '")\n'))
    cat(paste0("seqs_srcdir: ", path_to_multifasta, "\n"))
    cat(paste0("seqfiles_suffix: .fa.gz\n"))
    sink()
    c(package_name, seed_file_name)
}

# Get per-chromosome fasta file in .gz format:
generate_multifasta(path_single_fasta_genome = opt$genome_fasta, out_dir = "custom_genome_chromosome/")

# Get a seed file:
package_seed <- generate_seed_file(path_to_multifasta = "custom_genome_chromosome/", latin_name = opt$genome_name, common_name = "custom", genome_build = "custom", seed_file_name = "custom.BSgenome.seed.txt", fasta_url = "custo_url", version = "1.0.0")

# Forge the BSgenome package:
forgeBSgenomeDataPkg(package_seed[2], destdir = ".")

# Check, build and install:
dir.create("user_rlib", recursive = TRUE, showWarnings = FALSE)
.libPaths("user_rlib")  # add user rlib to the path

devtools::check(package_seed[1])
devtools::build(package_seed[1])
devtools::install(package_seed[1])

# Save to a .rds file:
library(BSgenome.Cgenome.Ensembl.custom)
saveRDS(BSgenome.Cgenome.Ensembl.custom, file = "custom_BSgenome.rds")
