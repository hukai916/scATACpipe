# Known issues

## 1. jpg image may not available in MultiQC html report
This might be due to too many genes/features are being plotted in a single jpg, which exceeds the jpg pixel limit.
If this happens, pls use the original pdf file.

## 2. Regarding pseudo-autosomal genes)  
When plotting, target genes are supplied with gene_symbols. for gene_ids that share the same gene symbol (pseudo-autosomal genes), all of them will be merged under the same gene symbol when plotting.

## 3. ArchR (Haibo repo) can't find parallel package
This might be due to package DESCRIPTION not configured properly. Fixed in Kai repo.


Dev notes:
For 1, one solution is to use png format, which supports much larger file, but the resulting images will be much larger making the html too heavy.
For 2, should always use gene_ids, however, ArchR is not implemented to take as input gene_ids when plotting.
For 3, a quick fix to the manually library(parallel).

For ARCHR_MOTIF_DEVIATIONS_CLUSTERS/CLUSTERS:
```
markerRNA <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
```
Above is okay for ArchR natively supported genome.
For Ensemble genome and other custom genomes, must exclude the leading "chrN:":
```
markerRNA <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- gsub(".+:", "", markerRNA) # get rid of potential leading "chrN:", otherwise plotEmbedding error
```
