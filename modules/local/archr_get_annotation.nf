// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Parse software version numbers
 */
process ARCHR_GET_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'archr_get_annotation', publish_id:'') }

    // conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    // } else {
    //     container "quay.io/biocontainers/python:3.8.3"
    // }

    // container "hukai916/bcl2fastq:2.20.0-centos7"
    container "hukai916/r_sc:0.5"

    // cache false

    input:
    val archr_genome

    output:
    path "genomeAnnotation.rds", emit: genomeAnnotation
    path "geneAnnotation.rds", emit: geneAnnotation
    path "user_rlib", emit: user_rlib

    script:

    """
    echo '
    library(ArchR)

    dict <- list(
      hg19 = c("BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"),
      hg38 = c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
      mm9 = c("BSgenome.Mmusculus.UCSC.mm9", "TxDb.Mmusculus.UCSC.mm9.knownGene", "org.Mm.eg.db"),
      mm10 = c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
      bosTau9 = c("BSgenome.Btaurus.UCSC.bosTau9", "TxDb.Btaurus.UCSC.bosTau9.refGene", "org.Bt.eg.db"),
      ce11 = c("BSgenome.Celegans.UCSC.ce11", "TxDb.Celegans.UCSC.ce11.refGene", "org.Ce.eg.db"),
      canFam3 = c("BSgenome.Cfamiliaris.UCSC.canFam3", "TxDb.Cfamiliaris.UCSC.canFam3.refGene", "org.Cf.eg.db"),
      danRer11 = c("BSgenome.Drerio.UCSC.danRer11", "TxDb.Drerio.UCSC.danRer11.refGene", "org.Dr.eg.db"),
      dm6 = c("BSgenome.Dmelanogaster.UCSC.dm6", "TxDb.Dmelanogaster.UCSC.dm6.ensGene", "org.Dm.eg.db"),
      galGal6 = c("BSgenome.Ggallus.UCSC.galGal6", "TxDb.Ggallus.UCSC.galGal6.refGene", "org.Gg.eg.db"),
      rheMac10 = c("BSgenome.Mmulatta.UCSC.rheMac10", "TxDb.Mmulatta.UCSC.rheMac10.refGene", "org.Mmu.eg.db"),
      panTro6 = c("BSgenome.Ptroglodytes.UCSC.panTro6", "TxDb.Ptroglodytes.UCSC.panTro6.refGene", "org.Pt.eg.db"),
      rn6 = c("BSgenome.Rnorvegicus.UCSC.rn6", "TxDb.Rnorvegicus.UCSC.rn6.refGene", "org.Rn.eg.db"),
      sacCer3 = c("BSgenome.Scerevisiae.UCSC.sacCer3", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "org.Sc.sgd.db"),
      susScr11 = c("BSgenome.Sscrofa.UCSC.susScr11", "TxDb.Sscrofa.UCSC.susScr11.refGene", "org.Ss.eg.db")
      )

    # Prepare genome/gene annotation
    if (!(length(dict["$archr_genome"]) == 3)) {
      dir.create("user_rlib", recursive = TRUE)
      .libPaths("user_rlib")  # add user rlib to the path

      if (!requireNamespace(dict["$archr_genome"][[1]][[1]], quietly = TRUE)) {
        BiocManager::install(dict["$archr_genome"][[1]][[1]])
      }

      if (!requireNamespace(dict["$archr_genome"][[1]][[2]], quietly = TRUE)){
        BiocManager::install(dict["$archr_genome"][[1]][[2]])
      }

      if (!requireNamespace(dict["$archr_genome"][[1]][[3]], quietly = TRUE)){
        BiocManager::install(dict["$archr_genome"][[1]][[3]])
      }

      library(dict["$archr_genome"][[1]][[1]], character.only = TRUE)
      library(dict["$archr_genome"][[1]][[2]], character.only = TRUE)
      library(dict["$archr_genome"][[1]][[3]], character.only = TRUE)

      genomeAnnotation <- createGenomeAnnotation(genome = get(dict["$archr_genome"][[1]][[1]]))
      geneAnnotation <- createGeneAnnotation(TxDb = get(dict["$archr_genome"][[1]][[2]]), OrgDb = get(dict["$archr_genome"][[1]][[3]]))

      saveRDS(genomeAnnotation, file = "genomeAnnotation.rds")
      saveRDS(geneAnnotation, file = "geneAnnotation.rds")

    } else {
      stop(paste0("Genome ", $archr_genome, " not supported!"))
    }

    ' > run.R

    Rscript run.R

    """
}
