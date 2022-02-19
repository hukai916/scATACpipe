# scATACpipe: Output

(Some parts adapted from nf-core [TEMPLATE](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/docs/output.md).)

## Table of Contents
[Introduction](#introduction)  
[Results](#results)  
  * [PREPROCESS_DEFAULT](#preprocess_default)
  * [PREPROCESS_10XGENOMICS](#preprocess_10xgenomics)
  * [PREPROCESS_CHROMAP](#preprocess_chromap)
  * [DOWNSTREAM_ARCHR](#downstream_archr)

[Examples](#example_reports)
  * [homo_sapiens](#)
  * [hg38](#)

[Pipeline information](#pipeline_information)

## Introduction
This document describes the output produced by scATACpipe. By default, all results are saved into `./results` folder unless otherwise specified by `--outdir` flag. Results from different Nextflow modules will be saved into corresponding folder (e.g. `./results/fastqc`, `./results/cutadapt`, *etc.*)

A html report summarizing all key results will be generated with MultiQC (via custom plugins) and saved into `./results/multiqc/multiqc_report.html` for quick view.

## Result folders
Nextflow implements a caching mechanism where all intermediate and final results are saved into `./work/` directory. By default, files under './results/' are symbolic links to that in the `./work/`. To switch from `symlink` to `copy`, use `--publish_dir_mode copy` flag

Below summarizes the main results by sub-workflow. The following modules are shared across multiple sub-workflows and are only described once under PREPROCESS_DEFAULT: `DOWNLOAD_FROM_UCSC`, `DOWNLOAD_FROM_UCSC_GTF`, `DOWNLOAD_FROM_ENSEMBL`, `DOWNLOAD_FROM_ENSEMBL_GTF`, `PREP_GENOME`, `PREP_GTF`.

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

#### PREPROCESS_CHROMAP

#### DOWNSTREAM_ARCHR

## Example reports
#### homo_sapiens
#### hg38

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
