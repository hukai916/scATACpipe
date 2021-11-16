library("GenomicFeatures")
# library("BSgenome.Hsapiens.Ensembl.GRCh38")
library("plyranges")
library("S4Vectors")
library("tibble")
library("dplyr")
library("biomaRt")
library("collections")
library("BSgenome")

# make TxDb from a GTF file
forgeTxDb <- function(BSgenome, gtf, out_TxDb_dir)
{
    if (missing(BSgenome) || missing(gtf) || missing(out_TxDb_dir))
    {
        stop("All arguments are required: BSgenome, gtf, and out_TxDb_dir!")
    }
    if (!is(BSgenome, "BSgenome")){
        stop("BSgenome is not a BSgenome object!")
    }
    if (!file.exists(gtf)){
        stop("gtf doesn't exist!")
    }
    if (!dir.exists(out_TxDb_dir)){
        dir.create(out_TxDb_dir, recursive = TRUE)
    }
    chrom_len <- seqlengths(BSgenome)
    is_circular <- names(chrom_len) %in% c("chrM", "chrMT", "MT",
                                           "chrPltd", "Pltd")
    chrominfo <- data.frame(chrom = names(chrom_len),
                            length = unname(chrom_len),
                            is_circular = is_circular)
    genome_metadata <- metadata(BSgenome)
    TxDb <- makeTxDbFromGFF(file = gtf,
                            format = "gtf",
                            dataSource = genome_metadata$provider,
                            organism = genome_metadata$organism,
                            taxonomyId = NA,
                            chrominfo = chrominfo,
                            miRBaseBuild = NA)
    TxDb_file <- file.path(out_TxDb_dir,
                           paste0(genome_metadata$genome, ".TxDb.sqlite"))
    saveDb(TxDb, file = TxDb_file)
    TxDb_file
}


## get gene Symbol from the GTF file if exists, otherwise get the gene symbols
## from biomart
get_geneID_symbol <- function(gtf = NULL, species_latin_name = NULL)
{
    if (is.null(gtf) || is.na(gtf) || !file.exists(gtf))
    {
        stop("gtf is required, or the gtf provided doesn't exist!")
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
    ## this can also be done using rtracklayer::import
    id2symbol_dict <- ordered_dict()
    gtf_full <- read.delim(in_gtf, header = FALSE, as.is = TRUE,
                      comment.char = "#", quote = "")
    close(in_gtf)
    gtf <- gtf_full[gtf_full[, 3] =="gene", 9]
    if (length(gtf) < 1)
    {
        message("There is not entries for genes in the gtf!\n", "Use transcript entries!")
        gtf <- gtf_full[gtf_full[, 3] =="transcript", 9]
        if (length(gtf) < 1) {stop("gtf is malformed! please doulbe check it!")}
    }
    null <- lapply(gtf, function(.x){
        if (grepl("gene_id", .x)){
            gene_id <- gsub('gene_id\\s+"([^".]+).+',
                            "\\1", .x, perl = TRUE)
            if (grepl("gene_name", .x))
            {
                gene_symbol <- gsub('.+?gene_name\\s+"([^"]+).*',
                                    "\\1", .x, perl = TRUE)
                if (grepl("^ENS.*?G\\d+", gene_symbol, perl = TRUE) ||
                    grepl("^ENS.*?T\\d+", gene_symbol, perl = TRUE) ||
                    grepl("^\\d+$", gene_symbol, perl = TRUE))
                {
                    gene_symbol <- "NA"
                }
            } else {
                gene_symbol <- "NA"
            }
            id2symbol_dict$set(gene_id, gene_symbol)
        }
    })

    ## check gene_id type
    gene_ids <- unlist(id2symbol_dict$keys())
    id_type <- {if (grepl("^ENS.*?G\\d+", gene_ids[1], perl = TRUE)) {"ensembl_gene_id"}
                else if (grepl("^\\d+$", gene_ids[1], perl = TRUE)){"entrezgene_id"}
                else if (!grepl("^ENS.*?T\\d+", gene_ids[1], perl = TRUE) &&
                         grepl("[a-zA-Z0-9]+", gene_ids[1],  perl = TRUE)) {"gene_name"}
                else {"unknown"}}
    if (id_type == "unknown")
    {
       stop("Unknown gene ID type!")
    } else if (id_type == "gene_name"){
       id_symbol <- data.frame(gene_id = gene_ids, symbol = gene_ids)
       return(id_symbol)
    }

    ## for gene_id types: ensembl_gene_id and entrezgene_id, if gene_name is sompletely missing, query Biomart databases
    if (all(unlist(id2symbol_dict$values()) == "NA"))
    {
        if (is.null(species_latin_name) || is.na(species_latin_name) || !grepl("^(.).+\\W(.+)", species_latin_name))
        {
            stop("species_latin_name is required or not correct!")
        }
        message("The GTF file contains no gene symbols. ",
                "Query Biomart to get gene symbols")
        require("biomaRt")
        species <- tolower(gsub("^(.).+\\W(.+)", "\\1\\2",
                              species_latin_name, perl = TRUE))

        ## try different biomart: animal, plant, fungi, metazona
        hosts <- c("https://www.ensembl.org/", "https://plants.ensembl.org/",
                   "https://fungi.ensembl.org/", "https://metazoa.ensembl.org/")
        marts <- try({marts <- sapply(hosts, function(.x){
                                 listMarts(host = .x)$biomart[1]})
                      marts
                   })
        while (is(marts, "try-error"))
        {
            marts <- try({marts <- sapply(hosts, function(.x){
                                     listMarts(host = .x)$biomart[1]})
                          marts
                        })
        }
        hosts <- marts
        id_symbol <- data.frame()
        for (i in seq_along(hosts))
        {
            datasets <- try({
                ensembl <- useEnsembl(biomart = hosts[i], host = names(hosts)[i])
                datasets <- searchDatasets(ensembl, pattern = species)$dataset
                datasets
                })
            while (is(datasets, "try-error"))
            {
                datasets <- try({
                                ensembl <- useEnsembl(biomart = hosts[i], host = names(hosts)[i])
                                datasets <- searchDatasets(ensembl, pattern = species)$dataset
                                datasets
                })
             }

            is_ds <- grepl(paste0(species, "_"), datasets)
            if (any(is_ds))
            {
                dataset <- datasets[is_ds]
                id_symbol <- try({
                    ensembl <- useMart(biomart = hosts[i],
                                   dataset = dataset,
                                   host = names(hosts)[i])
                   id_symbol <- biomaRt::select(ensembl, keys = gene_ids,
                                       columns = c(id_type,'external_gene_name'),
                                       keytype = id_type)
                   id_symbol
                })
                while (is(id_symbol, "try-error"))
                {
                    id_symbol <- try({
                        ensembl <- useMart(biomart = hosts[i],
                                           dataset = dataset,
                                           host = names(hosts)[i])
                        id_symbol <- biomaRt::select(ensembl, keys = gene_ids,
                                            columns = c(id_type,'external_gene_name'),
                                            keytype = id_type)
                        id_symbol
                     })
                  }

                # genes not in database
                unnamed_geneID_symbol <-
                    data.frame(id = gene_ids[!gene_ids %in% id_symbol[, 1]],
                               external_gene_name = "NA")
                colnames(unnamed_geneID_symbol)[1] <- id_type

                # no gene_name genes
                id_symbol$external_gene_name <- ifelse(id_symbol$external_gene_name =="", "NA", id_symbol$external_gene_name)
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
            stop("No gene symbols are get from Biomart!\n",
                    "Check your GTF file!")
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

## It is better to filter chrM, chrY for mammals, chrM for other animal and fungi, chrM and chrPltd for plants
create_ArchR_geneannotation_WO_OrgDb <- function(TxDb = NULL,
                                                 geneID2Symbol,
                                                 flank = 2000, ## prevent TSS expanded 2 Kb in both directions from being out of range
                                                 promoterRange = c(upstream = 2000, downstream = 2000), ## prevent extended promoters from being out of range; to be consistent with extended TSS, otherwise TSS lost while gene exists, which leads to error in create_ArchRgenome_Annotation.
                                                 filterChr = c("chrM", "chrY", "chrMT", "chrPltd"),
                                                 out_dir)
{
    if (is.null(TxDb) || !is(TxDb, "TxDb"))
    {
        stop("TxDb must ba a TxDb object!")
    }
    if (missing(geneID2Symbol) || missing(out_dir))
    {
        stop("geneID2Symbol and out are required!")
    }
    if (!is.data.frame(geneID2Symbol) ||
        any(colnames(geneID2Symbol) != c("gene_id", "symbol")))
    {
        stop("geneID2Symbol must be a dataframe with colnames:'gene_id', 'symbol'!")
    }
    if (!dir.exists(out_dir))
    {
        dir.create(out_dir, recursive = TRUE)
    }
    ## geneID2Symbol: col1: gene_id; col2: gene_symbol
    symbol <- geneID2Symbol[, 2]
    names(symbol) <- geneID2Symbol[, 1]

    ## filter seqnames of no interest, such as mitochondrial genome
    seqlevels_all <- seqlevels(TxDb)
    seqlevels(TxDb) <- seqlevels_all[!seqlevels_all %in% filterChr]

    ## GRanges for genes
    genes <- GenomicFeatures::genes(TxDb) %>%
        plyranges::remove_names()
    genes$symbol <- symbol[genes$gene_id]  # more stable
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

    ## get all transcripts
    tx <- unlist(transcriptsBy(TxDb, by = "gene")) %>%
        plyranges::mutate(gene = names(.)) %>%
        plyranges::remove_names() %>%
        plyranges::select(-c("tx_id")) %>%
        data.frame()
    tx_gene <- tx$gene
    names(tx_gene) <-tx$tx_name
    rm("tx")

    ## Create GRanges for exons
    exons <- unlist(exonsBy(TxDb,
                            by = "tx",
                            use.names = TRUE)) %>%
        plyranges::mutate(tx_name = names(.)) %>%
        plyranges::remove_names()
    exons$gene_id <- tx_gene[exons$tx_name]
    exons <- exons %>%
        plyranges::filter(!is.na(gene_id))
    exons$symbol <- symbol[exons$gene_id]
    exons <- exons %>%
        plyranges::select(-c("exon_id", "exon_name", "exon_rank",
                             "gene_id", "tx_name"))
    exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)

    ## Create GRanges for TSS
    TSS <- unique(resize(GenomicFeatures::transcripts(TxDb),
                         width = 1,
                         fix = "start")) %>%
           plyranges::select(-c("tx_id"))

    ## remove genes whose promoters are close to the chromosome end (promoter regions in upstream and downstream [2000, 100])
    gene_start <-  resize(genes, width = 1, fix ="start")
    gene_start_downstream <- stretch(anchor_5p(gene_start), extend = promoterRange[2])
    gene_start_upstream <- stretch(anchor_3p(gene_start), extend = promoterRange[1])
    downstream_out_of_bound_index <- GenomicRanges:::get_out_of_bound_index(gene_start_downstream)
    upstream_out_of_bound_index <- GenomicRanges:::get_out_of_bound_index(gene_start_upstream)
    gene_out_of_bound_index <- c(downstream_out_of_bound_index, upstream_out_of_bound_index)

    if (length(gene_out_of_bound_index) >0)
    {
        genes <- genes[-c(gene_out_of_bound_index)]
    }

    ## remove exons with genes removed due to out of bound
    exons <- exons[mcols(exons)$symbol %in% mcols(genes)$symbol]

    ## remove TSSs which are close to the chromosome end (<=2000 bp)
    TSS_2kb_flank <- resize(TSS,
                        width = 2 * flank + 1,
                        fix = "center")
    TSS_out_of_bound_index <- GenomicRanges:::get_out_of_bound_index(TSS_2kb_flank)
    if (length(TSS_out_of_bound_index) >0)  # otherwise get empty TSS
    {
        TSS <- TSS[-c(TSS_out_of_bound_index)]
    }

    ## drop unused seqlevels
    seqlevels(genes) <- seqlevelsInUse(genes)
    seqlevels(exons) <- seqlevelsInUse(exons)
    seqlevels(TSS) <- seqlevelsInUse(TSS)

    ## don't need ArchR createGeneAnnotation() function
    geneAnnotation <- SimpleList(genes = genes, exons = exons, TSS = TSS)
    saveRDS(geneAnnotation, file = file.path(out_dir, "geneAnnotation.RDS"))
    geneAnnotation
}

## create genomeAnnotation using a BSgenome, geneAnnotation, and an optional blacklist in the BED format.
## Always filter out chrM, chrPltd(plants), chrY (mammals), and chromosomes without TSSs/genes at all,
## Filter out chromosomes without TSSs/genes to avoid empty GRanges for TSS and genes after splitting
## by seqnames when create ArrowFiles which cause errors and no ArrowFiles will be generated for those excluded
## chromosomes. Usually, ArrowFiles for chromosomes without TSSs/genes are empty.
create_ArchR_genomeannotation <- function(BSgenome = NULL,
                                          geneAnnotation = NULL,
                                          out_dir = NULL,
                                          blacklist_bed = NULL,
                                          blacklist_hasheader = FALSE,
                                          filterChr = c("chrM", "chrY", "chrMT", "chrPltd"))
{
    if (missing(BSgenome) || is.null(BSgenome) || !is(BSgenome, "BSgenome"))
    {
        stop("BSgenome must be a BSgenomeobject!")
    }
    if (is.null(geneAnnotation) || !is(geneAnnotation, "SimpleList"))
    {
        stop("geneAnnotation must be a SimpleList containing GRanges for genes, exons and TSSs!")
    }
    if (is.null(out_dir))
    {
        stop("out_dir is required!")
    }
    if (!dir.exists(out_dir))
    {
        dir.create(out_dir, recursive = TRUE)
    }
    chrom_len <- seqlengths(BSgenome)
    chr_df <- data.frame(seqnames = names(chrom_len),
                         start = 1,
                         end = unname(chrom_len))
    chromSizes <- makeGRangesFromDataFrame(chr_df)
    all_seqlevels <- seqlevels(chromSizes)

    ## filter extra chromosomes/scaffolds so no ArrowFile generated for them
    tss_chr <- unique(as.character(seqnames(geneAnnotation$TSS)))
    if (!is.null(filterChr) && !is.na(filterChr))
    {
        tss_chr <- tss_chr[!tss_chr %in% filterChr]
    }
    ## filter chromSizes
    seqlevels(chromSizes, pruning.mode="coarse") <- tss_chr
    seqlevels(chromSizes) <- seqlevelsInUse(chromSizes)

    if (!is.null(blacklist_bed))
    {
        if (grepl(".bed.gz$", blacklist_bed))
        {
            blacklist <- gzfile(blacklist_bed, open ="rt")
        } else if (grepl(".bed$", blacklist_bed)) {
            blacklist <- file(blacklist_bed, open = "r")
        } else {
            stop("It seems the Blacklist file is not a bed file ",
                 "which should with an extension .bed, or .bed.gz")
        }
        blacklist_df <- read.delim(blacklist,
                                   header = blacklist_hasheader)
        close(blacklist)
        if (ncol(blacklist_df) < 3 ||
            any(!blacklist_df[,1] %in% names(chrom_len)) ||
            any(!is.numeric(blacklist_df[,2])) ||
            any(!is.numeric(blacklist_df[,3])))
        {
            stop("blacklist_bed is not a valid BED file!\n",
                 "Please make sure the chromosome names of the blacklist are a subset of those of the BSgenome.")
        }
        colnames(blacklist_df)[1:3] <- c("seqnames", "start", "end")
        blacklist <- makeGRangesFromDataFrame(blacklist_df,
                                              starts.in.df.are.0based = TRUE)
        ## filter blacklist
        seqlevels(blacklist, pruning.mode="coarse") <- tss_chr
        seqlevels(blacklist) <- seqlevelsInUse(blacklist)
    } else {
      blacklist <- GRanges()
    }

    # don't need ArchR createGenomeAnnotation() function
    genomeAnnotation <- SimpleList(genome = BSgenome@pkgname, chromSizes = chromSizes, blacklist = blacklist)
    saveRDS(genomeAnnotation, file = file.path(out_dir,
                                               "genomeAnnotation.RDS"))
    genomeAnnotation
}

# gene_symbol <-
#    get_geneID_symbol(gtf = "refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz",
#                      species_latin_name = "Homo sapiens")
