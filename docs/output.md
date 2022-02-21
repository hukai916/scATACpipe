# scATACpipe: Output

(Some parts adapted from nf-core [TEMPLATE](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/docs/output.md).)

## Table of Contents
[Introduction](#introduction)  
[Results](#results)  
  * [PREPROCESS_DEFAULT](#preprocess_default)
  * [PREPROCESS_10XGENOMICS](#preprocess_10xgenomics)
  * [PREPROCESS_CHROMAP](#preprocess_chromap)
  * [DOWNSTREAM_ARCHR](#downstream_archr)
  * [Other modules](#other_modules)

[Examples](#example_reports)
  * [Example1: default preprocessing with homo_sapiens (Ensembl)](#example1_default_preprocessing_with_homo_sapiens_ensembl))
  * [Example2: default preprocessing with hg38 (UCSC)](#)
  * [Example3: chromap preprocessing with homo_sapiens (Ensembl)](#)

[Pipeline information](#pipeline_information)

## Introduction
This document describes the output produced by scATACpipe. By default, all results are saved into `./results` folder unless otherwise specified by `--outdir` flag. Results from different Nextflow modules will be saved into corresponding folder (e.g. `./results/fastqc`, `./results/cutadapt`, *etc.*)

A html report summarizing all key results will be generated with MultiQC (via custom plugins) and saved into `./results/multiqc/multiqc_report.html` for quick view.

## Result folders
Nextflow implements a caching mechanism where all intermediate and final results are saved into `./work/` directory. By default, files under './results/' are symbolic links to that in the `./work/`. To switch from `symlink` to `copy`, use `--publish_dir_mode copy` flag

Below summarizes the main results for each sub-workflow. The following Nextflow modules are shared across multiple sub-workflows and are only described once under PREPROCESS_DEFAULT: `DOWNLOAD_FROM_UCSC`, `DOWNLOAD_FROM_UCSC_GTF`, `DOWNLOAD_FROM_ENSEMBL`, `DOWNLOAD_FROM_ENSEMBL_GTF`, `PREP_GENOME`, `PREP_GTF`.

Also see **module references**: [csv](https://github.com/hukai916/scATACpipe/blob/dev/docs/scATACpipe_module_references.csv) or [xlsx](https://github.com/hukai916/scATACpipe/blob/dev/docs/scATACpipe_module_references.xlsx) for more information.

#### PREPROCESS_DEFAULT
<details markdown="1">
<summary>Module: FASTQC</summary>

* `./results/fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).
</details>

<details markdown="1">
<summary>Module: SPLIT_FASTQ</summary>

* `./results/split_fastq/`
    * `*_xx.fastq.gz`: Chunked fastq files.
</details>

<details markdown="1">
<summary>Module: ADD_BARCODE_TO_READ(_CHUNKS)</summary>

* `./results/add_barcode_to_read(_chunks)/`
    * `*.barcoded.fastq.gz`: Barcode added fastq files.
* Prepend barcode to the front of the name line in FASTQ file. If `--split_fastq true`, output to `./results/add_barcode_to_read_chunks/`, else to `./results/add_barcode_to_read/`
</details>

<details markdown="1">
<summary>Module: CUTADAPT</summary>

* `./results/cutadapt/`
    * `*.trimmed.fastq.gz`: Adapter trimmed fastq files.
</details>

<details markdown="1">
<summary>Module: DOWNLOAD_FROM_ENSEMBL</summary>

* `./results/download_from_ensembl/`
    * `*.fa.gz`: Genome fasta downloaded from ENSEMBL.
</details>

<details markdown="1">
<summary>Module: DOWNLOAD_FROM_ENSEMBL_GTF</summary>

* `./results/download_from_ensembl_gtf/`
    * `*.gtf.gz`: GTF annotation downloaded from ENSEMBL.
</details>

<details markdown="1">
<summary>Module: DOWNLOAD_FROM_UCSC</summary>

* `./results/download_from_ucsc/`
    * `*.fa.gz`: Genome fasta downloaded from UCSC.
</details>

<details markdown="1">
<summary>Module: DOWNLOAD_FROM_UCSC_GTF</summary>

* `./results/download_from_ucsc_gtf/`
    * `*.gtf.gz`: GTF annotation downloaded from UCSC.
</details>

<details markdown="1">
<summary>Module: PREP_GENOME</summary>

* `./results/prep_genome/`
    * `*.fa.gz`: Cleaned genome fasta file.
* "Clean" means to extract primary genome and prepend "chr" in front of each chromosome.
</details>

<details markdown="1">
<summary>Module: PREP_GTF</summary>

* `./results/prep_gtf/`
    * `final.gtf.gz`: Cleaned GTF annotation file.
* "Clean" means to ensure that contigs in GTF are a subset of corresponding genome file; add prefix "chr" to each chromosome; and add "gene" column is not already there.
</details>

<details markdown="1">
<summary>Module: BWA_INDEX</summary>

* `./results/bwa_index/bwa_index`: BWA index files.
    * Genome.primary.chrPrefixed.fa.gz      
    * Genome.primary.chrPrefixed.fa.gz.bwt
    * Genome.primary.chrPrefixed.fa.gz.amb  
    * Genome.primary.chrPrefixed.fa.gz.pac
    * Genome.primary.chrPrefixed.fa.gz.ann  
    * Genome.primary.chrPrefixed.fa.gz.sa
</details>

<details markdown="1">
<summary>Module: BWA_MAP</summary>

* `./results/bwa_map/`: Mapped BAM files with BWA.
    * *.sorted.bam     
</details>

<details markdown="1">
<summary>Module: FILTER_BAM</summary>

* `./results/filter_bam/`: Filtered BAM files.
    * *.sorted.filtered.bam
* Filter out mitochondrial reads, and or 'improper' reads. To count as 'proper': paired reads must be mapped in the correct orientation; fragment size must range from 38 to 2000 bp; and the mapq of both reads must > 20.

</details>

<details markdown="1">
<summary>Module: COMBINE_BAM</summary>

* `./results/combine_bam/`: Combined BAM files.
    * *.sorted.bam
* Merge chunks/lanes of BAMs into a single BAM for each sample. This module will be used when `--barcode_correction false`.

</details>

<details markdown="1">
<summary>Module: COMBINE_BAM2</summary>

* `./results/combine_bam/`: Combined BAM files.
    * *.sorted.bam
* Merge chunks/lanes of BAMs into a single BAM for each sample. This module will be used when `--barcode_correction naive/pheniqs`.

</details>

<details markdown="1">
<summary>Module: DEDUP_BAM</summary>

* `./results/dedup_bam/`: Deduplicated BAM files using raw barcodes.
    * *.dedup.bam
* If `--barcode_correction false`, deduplicate BAM files using the raw barcodes.
</details>

<details markdown="1">
<summary>Module: DEDUP_BAM2</summary>

* `./results/dedup_bam/`: Deduplicated BAM files using barcode tags (corrected barcodes).
    * *.dedup.bam
* If `--barcode_correction naive/pheniqs`, deduplicate BAM files using the corrected barcodes.
</details>

<details markdown="1">
<summary>Module: GET_FRAGMENTS</summary>

* `./results/fragments/`: Extracted fragments from BAM files.
    * *.tsv.gz
</details>

<details markdown="1">
<summary>Module: QUALIMAP</summary>

* `./results/dedup/`: Various QC plots generated with Qualimap.
    * ./xxx/images_qualimapReport
    * ./xxx/qualimapReport.html
    * ./xxx/genome_results.txt
    * ./xxx/raw_data_qualimapReport
* `xxx` represents each sample name.
</details>

<details markdown="1">
<summary>Module: QUALIMAP</summary>

* `./results/dedup/`: Various QC plots generated with Qualimap.
    * ./xxx/images_qualimapReport
    * ./xxx/qualimapReport.html
    * ./xxx/genome_results.txt
    * ./xxx/raw_data_qualimapReport
* `xxx` represents each sample name.
</details>

#### PREPROCESS_10XGENOMICS
<details markdown="1">
<summary>Module: CELLRANGER_INDEX</summary>

* `./results/cellranger_index/genome_index`: Result folder for `cellranger-atac mkref`.
</details>

<details markdown="1">
<summary>Module: CELLRANGER_ATAC_COUNT</summary>

* `./results/cellranger_count`: Standard result folder for `cellranger-atac count`.
* More info regarding `cellranger-atac count` results, see [here](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count).
</details>

<details markdown="1">
<summary>Module: FILTER_CELL</summary>

* `./results/filter_cell/`
    * `*_valid_barcode_filtered_fragment.tsv.gz`: Filtered fragments by keeping only valid cells according to filtered_peak_bc_matrix/barcodes.tsv
    * `*_valid_barcode_filtered_fragment.bam`: Filtered BAM file by keeping only valid cells according to filtered_peak_bc_matrix/barcodes.tsv
</details>

#### PREPROCESS_CHROMAP
<details markdown="1">
<summary>Module: MERGE_SAMPLE</summary>

* `./results/merge_sample/`
    * `*.merge_read1.fastq.gz`: Lane-merged read1 fastq for each sample.
    * `*.merge_read2.fastq.gz`: Lane-merged read2 fastq for each sample.
    * `*.merge_barcode.fastq.gz`: Lane-merged barcode fastq for each sample.
</details>

<details markdown="1">
<summary>Module: GET_WHITELIST_CHROMAP</summary>

* `./results/get_whitelist_chromap/`
    * `*.txt.gz`: Selected whitelist barcodes.
</details>

<details markdown="1">
<summary>Module: CHROMAP_INDEX</summary>

* `./results/chromap_index/chromap_index_xxx`: Chromap genome index folder.
* `xxx` represents genome name.
</details>

<details markdown="1">
<summary>Module: CHROMAP_ATAC</summary>

* `./results/chromap_atac/*.sorted.tsv.gz`: Fragment file generated by `chromap --preset atac`.
</details>

<details markdown="1">
<summary>Module: FRAG_TO_FREQ</summary>

* `./results/frag_to_freq/`
    * `freq_xxx.sorted.tsv.gz`: Calculate fragment frequencies.
* `xxx` represents sample name.
</details>

<details markdown="1">
<summary>Module: GET_VALID_BARCODE_CHROMAP</summary>

* `./results/frag_to_freq/`
    * `xxx_valid_barcodes.txt`: Valid barcodes using "inflection point" method.
* `xxx` represents sample name.
</details>

<details markdown="1">
<summary>Module: FILTER_CELL_CHROMAP</summary>

* `./results/filter_cell_chromap/`
    * `xxx_valid_barcode_filtered_fragment.tsv.gz`: Valid fragments according to valid barcodes.
* `xxx` represents sample name.
</details>

#### DOWNSTREAM_ARCHR
<details markdown="1">
<summary>Module: PREP_FRAGMENT</summary>

* `./results/prep_fragment/*.fragment.tsv.gz`: Cleaned fragment file.
* `Clean` means to: prepend "chr" in front of fragment file to match PREP_GTF; ensure fragment file col1 is a subset of GTF col1.
</details>

<details markdown="1">
<summary>Module: ARCHR_CREATE_ARROWFILES</summary>

* `./results/archr_create_arrowfiles/`
    * `*.arrow`: Created arrowfiles using fragment files.
    * `QualityControl_`*: Quality Control results for created arrowfiles.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_CREATE_ARROWFILES_ANNOTATION</summary>

* `./results/archr_create_arrowfiles_annotation/`
    * `*.arrow`: Created arrowfiles using given custom annotation files.
    * `QualityControl_*`: Quality Control results for created arrowfiles.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_ADD_DOUBLETSCORES</summary>

* `./results/archr_add_doubletscores/`
    * `*._doublet.arrow`: Doublet score added arrowfiles.
    * `report_archr_add_doubletscores_pbmc_*`*: Reports on calculated doublet scores.
    * `summary_add_doubletscores_*.txt`
</details>

<details markdown="1">
<summary>Module: ARCHR_ARCHRPROJECT</summary>

* `./results/archr_archrproject/`
    * `ArchRProject`: Created ArchRProject folder.
    * `proj.rds`: A RDS file containing the ArchRProject object.
</details>

<details markdown="1">
<summary>Module: ARCHR_ARCHRPROJECT_ANNOTATION</summary>

* `./results/archr_archrproject_annotation/`
    * `ArchRProject`: Created ArchRProject folder using custom annotation files.
    * `proj.rds`: A RDS file containing the ArchRProject object.
</details>

<details markdown="1">
<summary>Module: ARCHR_ARCHRPROJECT_QC</summary>

* `./results/archr_archrproject_qc/`
    * `Plots`: QC plots.
    * `report_jpeg`: QC Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_GET_ANNOTATION_BIOC</summary>

* `./results/archr_get_annotation_bioc/`
    * `genomeAnnotation.rds`: Created ArchR genome annotation file.
    * `geneAnnotation.rds`: Created ArchR gene annotation file.
    * `user_rlib`: A directory containing installed custom R packages that can be passed to other modules.
</details>

<details markdown="1">
<summary>Module: BUILD_BSGENOME</summary>

* `./results/build_bsgenome/`
    * `custom_BSgenome.rds`: Created BSgenome.
    * `user_rlib`: A directory containing installed custom R packages that can be passed to other modules.
</details>

<details markdown="1">
<summary>Module: BUILD_TXDB</summary>

* `./results/build_txdb/`
    * `custom.TxDb.sqlite`: Created TxDb file.
</details>

<details markdown="1">
<summary>Module: BUILD_GENOME_ANNOTATION</summary>

* `./results/build_genome_annotation/`
    * `genomeAnnotation.RDS`: Created ArchR genome annotation file.
</details>

<details markdown="1">
<summary>Module: BUILD_GENE_ANNOTATION</summary>

* `./results/build_genome_annotation/`
    * `geneAnnotation.RDS`: Created ArchR gene annotation file.
</details>

<details markdown="1">
<summary>Module: ARCHR_FILTER_CELLS</summary>

* `./results/archr_filter_cells/`
    * `proj_cell_filtered.rds`: A RDS file containing cell-filtered ArchRProject.
    * Unqualified cells are filtered out by ArchR.
</details>

<details markdown="1">
<summary>Module: AMULET_DETECT_DOUBLETS</summary>

* `./results/amulet_detect_doublets/`
    * `cells_filter_*.txt`: AMULET detected doublet cells.
</details>

<details markdown="1">
<summary>Module: AMULET_MERGE_DOUBLETS</summary>

* `./results/amulet_merge_doublets/`
    * `amulet_doublets.txt`: Merged doublet cells for each sample.
</details>

<details markdown="1">
<summary>Module: AMULET_FILTER_DOUBLETS</summary>

* `./results/amulet_filter_doublets/`
    * `proj_doublet_filtered.rds`: A RDS file containing doublet filtered ArchRProj.
</details>

<details markdown="1">
<summary>Module: ARCHR_DIMENSION_REDUCTION</summary>

* `./results/archr_dimension_reduction/`
    * `proj_lsi.rds`: A RDS file containing ArchRProj with dimension reduction information.
</details>

<details markdown="1">
<summary>Module: ARCHR_BATCH_CORRECTION</summary>

* `./results/archr_batch_correction/`
    * `proj_batch_correct.rds`: A RDS file containing ArchRProj with batch correction information.
</details>

<details markdown="1">
<summary>Module: ARCHR_CLUSTERING</summary>

* `./results/archr_clustering/`
    * `proj_clustering.rds`: A RDS file containing ArchRProj with clustering information.
    * `Clusters_matrix.csv`: A CSV file containing clustering information.
    * `Plots`: Clustering plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_EMBEDDING</summary>

* `./results/archr_embedding/`
    * `proj_embedding.rds`: A RDS file containing ArchRProj with embedding information.
    * `Plots`: Embedding plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_PSEUDO_BULK_CLUSTERS</summary>

* `./results/archr_pseudo_bulk_clusters/`
    * `archr_project_pseudobulk.rds`: A RDS file containing ArchRProj with pseudobulk information.
    * `archr_project_pseudobulk`: A ArchRproject folder with pseudobulk information.
    * `user_rlib`: A directory containing installed custom R packages that can be passed to other modules.
</details>

<details markdown="1">
<summary>Module: ARCHR_PSEUDO_BULK_CLUSTERS2</summary>

* `./results/archr_pseudo_bulk_clusters2/`
    * `archr_project_pseudobulk.rds`: A RDS file containing ArchRProj with pseudobulk information.
    * `archr_project_pseudobulk`: A ArchRproject folder with pseudobulk information.
    * `user_rlib`: A directory containing installed custom R packages that can be passed to other modules.
</details>

<details markdown="1">
<summary>Module: ARCHR_SCRNASEQ_UNCONSTRAINED</summary>

* `./results/archr_scrnaseq_unconstrained/`
    * `proj_scrnaseq_unconstrained.rds`: A RDS file containing ArchRProj with unconstrained scRNAseq integration information.
    * `cell_type_scRNA.txt`: cell type information stored inside scRNAseq data.
    * `Plots`: Unconstrained scRNAseq integration plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_SCRNASEQ_CONSTRAINED</summary>

* `./results/archr_scrnaseq_constrained/`
    * `proj_scrnaseq_constrained.rds`: A RDS file containing ArchRProj with unconstrained scRNAseq integration information.
    * `Plots`: Constrained scRNAseq integration plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_MARKER_GENE_CLUSTERS</summary>

* `./results/archr_marker_gene_clusters/`
    * `proj_marker_gene.rds`: A RDS file containing ArchRProj with marker gene information.
    * `Plots`: Marker gene plots.
    * `report_jpeg`: Plots in jpeg format.
    * `marker_list.txt`: Marker gene list.
    * `Clusters_markerList.rds`: A RDS file containing markerList information.
</details>

<details markdown="1">
<summary>Module: ARCHR_MARKER_GENE_CLUSTERS2</summary>

* `./results/archr_marker_gene_clusters2/`
    * `proj_marker_gene.rds`: A RDS file containing ArchRProj with marker gene information.
    * `Plots`: Marker gene plots.
    * `report_jpeg`: Plots in jpeg format.
    * `marker_list.txt`: Marker gene list.
    * `Clusters2_markerList.rds`: A RDS file containing markerList information.
</details>

<details markdown="1">
<summary>Module: ARCHR_CALL_PEAKS_CLUSTERS</summary>

* `./results/archr_call_peaks_clusters/`
    * `proj_call_peaks.rds`: A RDS file containing ArchRProj with called peaks.
</details>

<details markdown="1">
<summary>Module: ARCHR_CALL_PEAKS_CLUSTERS2</summary>

* `./results/archr_call_peaks_clusters2/`
    * `proj_call_peaks.rds`: A RDS file containing ArchRProj with called peaks.
</details>

<details markdown="1">
<summary>Module: ARCHR_GET_MARKER_PEAKS_CLUSTERS</summary>

* `./results/archr_get_marker_peaks_clusters/`
    * `proj_call_peaks.rds`: A RDS file containing ArchRProj with called peaks.
    * `Plots`: Marker peak plots.
    * `report_jpeg`: Plots in jpeg format.
    * `Clusters_group_names.txt`: Cluster group names.
    * `Clusters_marker_peaks.rds`: A RDS file containing marker peaks information.
</details>

<details markdown="1">
<summary>Module: ARCHR_GET_MARKER_PEAKS_CLUSTERS2</summary>

* `./results/archr_get_marker_peaks_clusters2/`
    * `proj_call_peaks.rds`: A RDS file containing ArchRProj with called peaks.
    * `Plots`: Marker peak plots.
    * `report_jpeg`: Plots in jpeg format.
    * `Clusters2_group_names.txt`: Cluster group names.
    * `Clusters2_marker_peaks.rds`: A RDS file containing marker peaks information.
</details>

<details markdown="1">
<summary>Module: ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS</summary>

* `./results/archr_marker_peaks_in_tracks_clusters/`
    * `Plots`: Track plots for marker peaks.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_MARKER_PEAKS_IN_TRACKS_CLUSTERS2</summary>

* `./results/archr_marker_peaks_in_tracks_clusters2/`
    * `Plots`: Track plots for marker peaks.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_PAIRWISE_TEST_CLUSTERS</summary>

* `./results/archr_pairwise_test_clusters/`
    * `test_group.txt`: Which groups to test against.
    *  `markerTest.rds`: A RDS file containing getMarkerFeatures() results.
    * `Plots`: Track plots for marker peaks.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_PAIRWISE_TEST_CLUSTERS2</summary>

* `./results/archr_pairwise_test_clusters2/`
    * `test_group.txt`: Which groups to test against.
    *  `markerTest.rds`: A RDS file containing getMarkerFeatures() results.
    * `Plots`: Track plots for marker peaks.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_MOTIF_ENRICHMENT_CLUSTERS</summary>

* `./results/archr_motif_enrichment_clusters/`
    * `archr_motif_enrichment_project.rds`: A RDS file containing motif enrichment information.
    * `Plots`: Motif enrichment plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_MOTIF_ENRICHMENT_CLUSTERS2</summary>

* `./results/archr_motif_enrichment_clusters2/`
    * `archr_motif_enrichment_project.rds`: A RDS file containing motif enrichment information.
    * `Plots`: Motif enrichment plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_MOTIF_DEVIATIONS_CLUSTERS</summary>

* `./results/archr_motif_deviations_clusters/`
    * `archr_motif_deviation_project.rds`: A RDS file containing motif deviation information.
    * `motif_names.txt`
    * `Plots`: Motif deviation plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_MOTIF_DEVIATIONS_CLUSTERS2</summary>

* `./results/archr_motif_deviations_clusters2/`
    * `archr_motif_deviation_project.rds`: A RDS file containing motif deviation information.
    *  `motif_names.txt`
    * `Plots`: Motif deviation plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_FOOTPRINTING_CLUSTERS</summary>

* `./results/archr_footprinting_clusters/`
    * `save_archr_project`: A ArchRproject folder with footprinting information.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_FOOTPRINTING_CLUSTERS2</summary>

* `./results/archr_footprinting_clusters2/`
    * `save_archr_project`: A ArchRproject folder with footprinting information.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_COACCESSIBILITY_CLUSTERS</summary>

* `./results/archr_coaccessibility_clusters/`
    * `archr_coaccessibility_project.rds`: A RDS file containing coaccessibility information.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_COACCESSIBILITY_CLUSTERS2</summary>

* `./results/archr_coaccessibility_clusters2/`
    * `archr_coaccessibility_project.rds`: A RDS file containing coaccessibility information.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_PEAK2GENELINKAGE_CLUSTERS2</summary>

* `./results/archr_coaccessibility_clusters2/`
    * `Plots`: Coaccessibility plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS</summary>

* `./results/archr_get_positive_tf_regulator_clusters/`
    * `Plots`: Positive TF regulator plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_GET_POSITIVE_TF_REGULATOR_CLUSTERS2</summary>

* `./results/archr_get_positive_tf_regulator_clusters2/`
    * `Plots`: Positive TF regulator plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_TRAJECTORY_CLUSTERS2</summary>

* `./results/archr_trajectory_clusters2/`
    * `Plots`: Trajectory plots.
    * `report_jpeg`: Plots in jpeg format.
</details>

<details markdown="1">
<summary>Module: ARCHR_GET_CLUSTERING_TSV</summary>

* `./results/archr_get_clustering_tsv/`
    * `Clusters*.tsv`: A TSV file containing cell clustering inforamtion.
</details>

#### Other modules
<details markdown="1">
<summary>Module: SPLIT_BED</summary>

* `./results/split_bed/`
    * `split_Clusters_*/cluster_*.txt`: Fragments by cluster.
</details>

<details markdown="1">
<summary>Module: SPLIT_BAM</summary>

* `./results/split_bam/`
    * `split_Clusters*/*.bam`: BAM files by cluster.
    * `bigWig_Clusters*/*.bw`: bigWig files by cluster.
</details>

<details markdown="1">
<summary>Module: MULTIQC</summary>

* `./results/multiqc/`
    * `*`: Cached result files/folders as input to MultiQC.
    * `multiqc_report.html`: A MultiQC report integration all core results.
</details>

## Examples

#### Example1: default preprocessing with homo_sapiens (Ensembl)
Command used:
```bash
cd scATACpipe
nextflow run main.nf --outdir test_data1_default_homo_sapiens --input_fastq assets/sample_sheet_test_data1.csv --preprocess default --ref_fasta_ensembl homo_sapiens --species_latin_name 'homo sapiens' --doublet_removal_algorithm amulet --amulet_rmsk_bed assets/rmsk_hg38.bed  --amulet_autosomes assets/homo_sapiens_autosomes.txt -profile singularity,lsf
```
Specifically:  
```
--outdir test_data1_default_homo_sapiens
```
Above means that all results to be saved into `test_data1_default_homo_sapiens directory`.

```
--ref_fasta_ensembl homo_sapiens
```
Above means that the Ensembl genome `homo_sapiens` will be downloaded and served as reference.

```
--preprocess default
```
Above means that the default preprocessing strategy will be used.

```
--doublet_removal_algorithm amulet
```
Above means that AMULET will be used for doublet removal.

```
-profile singularity,lsf
```
Above means Singularity will be used and jobs will be submitted via LSF.

A final report is saved at [./test_data1_default_homo_sapiens/multiqc/multiqc_report.html](#https://rawcdn.githack.com/hukai916/scATACpipe_example/2f396d2abcfcedd902bee4e6b9cb8674e2181d57/test_data1_default_homo_sapiens/multiqc/multiqc_report.html)

#### Example2: default preprocessing with hg38 (UCSC)
Command used:
```bash
cd scATACpipe
nextflow run main.nf --outdir test_data1_default_hg38 --input_fastq assets/sample_sheet_test_data1.csv --preprocess default --ref_fasta_ucsc hg38 -profile singularity,lsf
```
Specifically:  
```
--outdir test_data1_default_hg38
```
Above means that all results to be saved into `test_data1_default_hg38 directory`.

```
--ref_fasta_ucsc hg38
```
Above means that the UCSC genome `hg38` will be downloaded and served as reference.

```
--preprocess default
```
Above means that the default preprocessing strategy will be used. Since `--doublet_removal_algorithm` is not specified, the default `--doublet_removal_algorithm archr` will be used.

```
-profile singularity,lsf
```
Above means Singularity will be used and jobs will be submitted via LSF.

A final report is saved at [./test_data1_default_hg38/multiqc/multiqc_report.html](#https://rawcdn.githack.com/hukai916/scATACpipe_example/2f396d2abcfcedd902bee4e6b9cb8674e2181d57/test_data1_default_hg38/multiqc/multiqc_report.html)

#### Example3: preprocessing with Chromap using hg38 (UCSC)
Command used:
```bash
cd scATACpipe
nextflow run main.nf --outdir test_data1_chromap_hg38 --input_fastq assets/sample_sheet_test_data1.csv --preprocess chromap --barcode_correction naive --ref_fasta_ucsc hg38 --species_latin_name 'homo sapiens' -profile singularity,lsf
```
Specifically:  
```
--outdir test_data1_chromap_hg38
```
Above means that all results to be saved into `test_data1_chromap_hg38 directory`.

```
--ref_fasta_ucsc hg38
```
Above means that the UCSC genome `hg38` will be downloaded and served as reference.

```
--preprocess chromap
```
Above means that the Chromap preprocessing strategy will be used.

```
--barcode_correction naive
```
Above means that the naive barcode correction strategy will be used.

```
-profile singularity,lsf
```
Above means Singularity will be used and jobs will be submitted via LSF.

A final report is saved at [./test_data1_chromap_hg38/multiqc/multiqc_report.html](#https://rawcdn.githack.com/hukai916/scATACpipe_example/2f396d2abcfcedd902bee4e6b9cb8674e2181d57/test_data1_chromap_hg38/multiqc/multiqc_report.html)

## Pipeline information

A MultiQC report (.html) summarizing all key results

will be generated by


When inputting FASTQ files (`--input_fastq`), scATACpipe provides 3 alternative strategies for preprocessing (shown below), followed by DOWNSTREAM_ARCHR.
  - PREPROCESS_DEFAULT
  - PREPROCESS_CHROMAP
  - PREPROCESS_10XGENOMICS

When inputting Fragment files (`--input_fragment`), DOWNSTREAM_ARCHR will be directly executed.

The 2 input options are demonstrated in details below. You can also view all available parameters by following [Quick Start](https://github.com/hukai916/scATACpipe/#quick-start) and run:
```bash
cd scATACpipe
nextflow run main.nf --help
```

## Introduction

This document describes the output produced by scATACpipe. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [FastQC](#fastqc) - Raw read QC
* [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.tsv`.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
