#!/usr/bin/env Rscript

library(optparse)
option_list = list(
	make_option(c("--gtf"), type="character", default=NULL,
							help="full path to gtf file", metavar="character"),
	make_option(c("--species_latin_name"), type="character", default=NULL,
							help="Species latin name (e.g. 'Homo Sapiens')", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$gtf)) {
	print_help(opt_parser)
	stop("The gtf (.gtf.gz or .gtf) must be provided!")
}

library(dplyr)
library(biomaRt)
library(collections)

## get gene Symbol from the GTF file if exists, otherwise get the gene symbols
## from biomart
get_geneID_symbol <- function(gtf, species_latin_name)
{
    if (!file.exists(gtf)){
        stop(gtf, " doesn't exist!")
    }
    if (grepl(".gtf.gz$", gtf))
    {
        in_gtf <- gzfile(gtf, open ="r")
    } else if (grepl(".gtf$", gtf)) {
        in_gtf <- file(gtf, open = "r")
    } else {
        stop("It seems the GTF file is not a GTF file ",
             "which should with an extension .gtf, or .gtf.gz")
    }
    require("collections")
    id2symbol_dict <- ordered_dict()
    gtf <- read.delim(in_gtf, header = FALSE, as.is = TRUE,
                      comment.char = "#", quote = "")
    close(in_gtf)
    gtf <- gtf[gtf[, 3] == "gene", 9]
    null <- lapply(gtf, function(.x){
        if (grepl("gene_id", .x)){
            gene_id <- gsub('gene_id\\s+"([^"]+).+',
                            "\\1", .x, perl = TRUE)
            if (grepl("gene_name", .x))
            {
                gene_symbol <- gsub('.+?gene_name\\s+"([^"]+).+',
                                    "\\1", .x, perl = TRUE)
            } else {
                gene_symbol <- "NA"
            }
            id2symbol_dict$set(gene_id, gene_symbol)
        }
    })

    if (all(unlist(id2symbol_dict$values()) == "NA"))
    {
        message("The GTF file contains no gene symbols. ",
                "Query Biomart to get gene symbols")
        require("biomaRt")
        species <- tolower(gsub("^(.).+\\W(.+)", "\\1\\2",
                              species_latin_name, perl = TRUE))

        ## try different biomart: animal, plant, fungi, metazona
        hosts <-gene_ids <- unlist(id2symbol_dict$keys())
        id_type <- {if (grepl("^ENS.*?G\\d+", gene_ids[1])) {"ensembl_gene_id"}
            else if (grepl("^\\d+$", gene_ids[1])){"entrezgene_id"}
            else {"unknown"}}
        if (id_type == "unknow")
        {
            stop("Unkow gene ID type!")
        }


        hosts <- c("https://www.ensembl.org/", "https://plants.ensembl.org/",
                   "https://fungi.ensembl.org/", "https://metazoa.ensembl.org/")

        hosts <- sapply(hosts, function(.x){
            listMarts(host = .x)$biomart[1]
        })

        id_symbol <- data.frame()
        for (i in seq_along(hosts))
        {
            ensembl <- useEnsembl(biomart = hosts[i], host = names(hosts)[i])
            datasets <- searchDatasets(ensembl, pattern = species)$dataset
            is_ds <- grepl(paste0(species, "_"), datasets)
            if (any(is_ds))
            {
                dataset <- datasets[is_ds]
                ensembl <- useMart(biomart = hosts[i],
                                   dataset = dataset,
                                   host = names(hosts)[i])
                id_symbol <- select(ensembl, keys = gene_ids,
                       columns = c(id_type,'external_gene_name'),
                       keytype = id_type)

                unnamed_geneID_symbol <-
                    data.frame(id = gene_ids[!gene_ids %in% id_symbol[ ,1]],
                               external_gene_name = "NA")
                colnames(unnamed_geneID_symbol)[1] <- id_type
                id_symbol <- rbind(id_symbol, unnamed_geneID_symbol)
                id_symbol$external_gene_name <- paste(id_symbol[, 2],
                                                      id_symbol[, 1],
                                                      sep = "_")
                break
            }
        }
        if (nrow(id_symbol) > 0){
            colnames(id_symbol) <- c("gene_id", "symbol")
            return(id_symbol)
        } else {
            message("No gene symbols are get from Biomart!\n",
                    "Check your GTF file!")
            return(NULL)
        }
    } else {
        id2symbol <- data.frame(gene_id = unlist(id2symbol_dict$keys()),
                                     symbol = unlist(id2symbol_dict$values()))
        id2symbol$symbol <- paste(id2symbol$symbol,
                                       id2symbol$gene_id,
                                       sep = "_")
        return(id2symbol)
    }
}

gtf <- opt$gtf
species_latin_name <- opt$species_latin_name

id_symbol <- get_geneID_symbol(gtf = gtf, species_latin_name = species_latin_name)
saveRDS(id_symbol, file = "id_symbol.rds")
