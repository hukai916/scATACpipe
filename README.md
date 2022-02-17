# ![scATACpipe](docs/images/scATACpipe.png)

## Table of Contents
[Introduction](#introduction)

[Pipeline summary](#pipeline-summary)
- [PREPROCESS_DEFAULT](#preprocess_default)
- [PREPROCESS_10XGENOMICS](#preprocess_10xgenomics)
- [PREPROCESS_CHROMAP](#preprocess_chromap)
- [DOWNSTREAM_ARCHR](#downstream_archr)

[Quick Start](#quick-start)

[Documentation](#documentation)

[Credits](#credits)

[Bug report/Support](#bug_report)

[Citations](#citations)

## Introduction

**scATACpipe** is a bioinformatic pipeline for single-cell ATAC-seq (scATAC-seq) data analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity containers making installation trivial and results highly reproducible.

The development of the pipeline is guided by  [nf-core TEMPLATE](https://github.com/nf-core/tools/tree/master/nf_core/pipeline-template).

## Pipeline Summary

The pipeline consists of 2 relevant parts: preprocessing (from fastq to fragment file) and downstream analysis. If fragment files are directly available, you can choose to skip preprocessing and run downstream analysis only.

For preprocessing, 3 alternative strategies are available that are implemented in 3 sub-workflows respectively, namely, **PREPROCESS_DEFAULT**, **PREPROCESS_10XGENOMICS**, and **PREPROCESS_CHROMAP**. Each of them supports various input types that are demonstrated in further detail below (also see [usage](https://github.com/hukai916/scATACpipe/blob/main/docs/usage.md)).

For downstream analysis, we implemented **DOWNSTREAM_ARCHR** sub-workflow that integrates ArchR and other tools (e.g. AMULET for doublet detection).

Below is a simplified diagram to illustrate the design logic and functionalities of scATACpipe.

![scATACpipe](docs/images/scATACpipe_workflow.svg)

The main functionalities of each sub-workflow are summarized below:

#### PREPROCESS_DEFAULT:
1. Add barcodes to reads
2. Correct barcodes (optional)
    * meanwhile, also filter out non-cells
3. Trim off adapters
4. Mapping
    * download genome/annotation or use custom genome
    * build genome index if not supplied
5. Filter BAM
6. Remove PCR duplicates
7. Quality control
8. Generate fragment file, *etc.*

#### PREPROCESS_10XGENOMICS:
1. Build 10XGENOMICS index if not supplied
    * download genome/annotation or use custom genome
2. Execute `cellranger_atac count` command

#### PREPROCESS_CHROMAP:
1. Build Chromap index if not supplied
    * download genome/annotation or use custom genome
2. Execute `chromap --preset atac` command
3. Filter out non-cells

Note that no BAM file will be generated for PREPROCESS_CHROMAP option.

#### DOWNSTREAM_ARCHR:
1. Build ArchR-compatible genome/annotation files if not natively supported (ArchR supports hg19, hg38, mm9, and mm10 as of 02/2022)
    * download genome/annotation if not supplied
    * build ArchR genome/gene annotation files if needed
2. Perform downstream analysis with ArchR and generate various analytical plots
    * filter doublets (with ArchR built-in method or AMULET)
    * dimension reduction
    * batch effect correction
    * clustering
    * embedding
    * pseudo-bulk clustering
    * scRNAseq integration if supplied
    * marker gene detection
    * call peaks
    * marker peak detection
    * pairwise testing
    * motif enrichment
    * footprinting
    * coaccessibility, *etc.*

The pipeline also splits BED and/or BAM files according to ArchR clusterings and summarizes all results into a single MultiQC report for easy view.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)(>=21.10.0).

2. For full reproducibility, install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity (Apptainer)`](https://www.sylabs.io/guides/3.0/user-guide/) and specify `-profile singularity` or `-profile docker` accordingly when running the pipeline so that all dependencies are satisfied. Otherwise, all of the dependencies must be available locally on your PATH, **which is likely not true!**

3. Download the pipeline:
```bash
git clone https://github.com/hukai916/scATACpipe.git
```

4. Download a minimal test dataset:
    * The **test_data1** is prepared by downsampling (5% and 10%) a dataset named "*500 Peripheral blood mononuclear cells (PBMCs) from a healthy donor (Next GEM v1.1)*" provided by [10xgenomics](https://www.10xgenomics.com/resources/datasets?query=&page=1&configure%5Bfacets%5D%5B0%5D=chemistryVersionAndThroughput&configure%5Bfacets%5D%5B1%5D=pipeline.version&configure%5BhitsPerPage%5D=500&menu%5Bproducts.name%5D=Single%20Cell%20ATAC). Note that, in test_data1, I1 refers to index1, which is for sample demultiplexing and not relevant in our case; R1 refers to Read1; **R2 refers to index2**, which represents the cell barcode fastq; R3 refers to Read2.

  ```bash
  cd scATACpipe
  wget https://www.dropbox.com/s/uyiq18zk7dts9fx/test_data1.zip
  unzip test_data1.zip
  ```


5. Edit the `replace_with_full_path` in the assets/sample_sheet_test_data1.csv to use the actual **full path**.

6. Test the pipeline with this minimal test_data1:
    * At least 8GB memory is recommended for test_data1.
    * By default, the local executor will be used (`-profile local`) meaning that all jobs will be executed on your local computer.  Nextflow supports many other [executors](https://www.nextflow.io/docs/latest/executor.html) including SLURM, LSF, *etc.*. You can create a [profile](https://www.nextflow.io/docs/latest/config.html?highlight=profile#config-profiles) file to config which executor to use. Multiple profiles can be supplied with comma, e.g. `-profile docker,lsf`.
    * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see what other custom config files can be supplied.

  * **Example run with Docker:**
  ```bash
  nextflow run main.nf -profile docker --preprocess default --outdir res_test_data1 --input_fastq assets/sample_sheet_test_data1.csv --ref_fasta_ensembl homo_sapiens --species_latin_name 'homo sapiens'
  ```
    By executing the above command:
      - The `local executor` will be used.
      - `PREPROCESS_DEFAULT` will be used.
      - Output will be saved into `res_test_data1`.
      - Ensembl genome `homo_sapiens` will be downloaded and used as reference.

 * **Example run with Singularity:**
  ```bash
  nextflow run main.nf -profile singularity,lsf --preprocess default --outdir res_test_data1 --input_fastq assets/sample_sheet_test_data1.csv --ref_fasta_ensembl homo_sapiens --species_latin_name 'homo sapiens'
  ```
      - By specifying `-profile lsf`, the `lsf` executor will be used for job submission.
      - By specifying `-profile singularity`, Singularity images will be downloaded and saved to `work/singularity` directory. It is recommended to config the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) settings to store the images in a central location.

7. Run your own analysis:
  * A typical command:
```bash
nextflow run main.nf -profile <singularity/docker/lsf> --preprocess <default/10xgenomics/chromap> --outdir <path_to_result_dir> --input_fastq <path_to_samplesheet> --ref_fasta_ensembl <ENSEMBL_genome_name> --species_latin_name <e.g. 'homo sapiens'>
```
  * For help:
```bash
nextflow run main.nf --help
```

See documentation [usage](https://github.com/hukai916/scATACpipe/blob/main/docs/usage.md) for all of the available options.

## Documentation

The scATACpipe workflow comes with documentation about the pipeline: [usage](https://github.com/hukai916/scATACpipe/blob/main/docs/usage.md) and [output](https://github.com/hukai916/scATACpipe/blob/main/docs/output.md).

## Credits

scATACpipe was originally designed and written by Kai Hu, Haibo Liu, and Julie Lihua Zhu.

We thank the following people for their extensive assistance in the development
of this pipeline: Nathan Lawson.

## Bug report/Support

For help, bug report, or feature requests, the developers would prefer and appreciate that you create a GitHub issue by clicking [here](https://github.com/hukai916/scATACpipe/issues/new/choose).
If you would like to extend scATACpipe for your own good, feel free to fork the repo.

## Citations
<!-- TODO If you use scATACpipe for your analysis, please cite it using the following doi: [](https://) -->
Please kindly cite scATACpipe [to be added] if you use it for your research.

A complete list of references for the tools used by scATACpipe can be found [here](https://github.com/hukai916/scATACpipe/docs/module_software.xlsx). Please consider citing them too, it won't cost a penny.
