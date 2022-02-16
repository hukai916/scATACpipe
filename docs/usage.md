# scATACpipe: Usage

## Introduction
Either raw sequencing data (**.fastq**) or preprocessed fragment (**.gz**) file can serve as input to scATACpipe.

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

For a complete list of implemented scATACpipe modules, see **module references**: [csv](https://github.com/hukai916/scATACpipe/blob/dev/docs/scATACpipe_module_references.csv) or [xlsx](https://github.com/hukai916/scATACpipe/blob/dev/docs/scATACpipe_module_references.xlsx).

### Basics:
```
--input_fragment        [string]  Path to input sample sheet for fragment files.
--input_fastq           [string]  Path to input sample sheet for FASTQ files.
--outdir                [string]  Path to result folder. Default to ./results.
--support_genome        Show currently supported genomes.
```
To view currently support genomes, simply:
```bash
cd scATACpipe
nextflow run main.nf --support_genome
```
Refer to [output]() for example commands and results.

## Fragment files as input
Fragment file paths (full path) must be saved into a **.csv** file (see below) and supplied with `--input_fragment`.
```
sample_name,file_path            
SAMPLE_1,/full_path/xxx.tsv.gz
SAMPLE_2,/full_path/xxx.tsv.gz
```
| Column       | Description                                                |
|--------------|------------------------------------------------------------|
| sample_name  | Sample name must not contain space in it.                  |
| file_path    | Must use full path. File must have the extension '.tsv.gz'.|

An example .csv can be found [here](https://raw.githubusercontent.com/hukai916/scATACpipe/main/assets/example_samplesheet_fragment.csv).

In addition to the fragment files, genome/annotation files must also be supplied and there are 3 options, see below.

### Option 1: using UCSC/ENSEMBL genome
```
--archr_genome          [string]  A genome name, either ENSEMBL style (e.g. homo_sapiens) or UCSC style (e.g. mm10).
--species_latin_name    [string]  Must be quoted. Required if '--archr_genome' not in (mm9, mm10, hg19, hg38)
--archr_blacklist       [string]  Optional. Path to blacklist file.
```

### Option2: using custom genome
```
--archr_genome_fasta    [string]  Path to genome fasta.
--ref_gtf               [string]  Path to gtf file.
--species_latin_name    [string]  Must be quoted.
--archr_blacklist       [string]  Optional. Path to blacklist file.
```

### Option3: using Bioconductor annotations
```
--archr_bsgenome    [string]  A Bioconductor BSgenome package name.
--archr_txdb        [string]  A Bioconductor TxDb package name.
--archr_org         [string]  A Bioconductor OrgDb package name.
--archr_blacklis    [string]  Optional. Path to blacklist file.
```

### Other parameters
Given that downstream analysis itself is highly interactive in nature, scATACpipe was implemented in a way that is as flexible as possible, meaning that users can configure many downstream parameters.

The parameters can be divided into two categories, namely, **main pipeline parameters**, and **module specific parameters**.

Main pipeline parameters must be supplied with command flags or configured inside `nextflow.confg`. They are typically required to instruct scATACpipe to perform certain analysis. These parameters are listed below.
```
--archr_thread                      [Integer]  Number of threads to use. Default to 4.

--archr_batch_correction_harmony    [true|false]  Whether or not to perform batch correction with Harmony.

--doublet_removal_algorithm         [false|archr|amulet]  Doublet removal algorithm, use false to skip.
  # If "amulet", also need to config: conf/modules.config -> amulet_detect_doublets.args and the followings:
--amulet_rmsk_bed                   [string]  Path to known repeat regions in genome.
--amulet_autosomes                  [string]  Path to txt file containing autosomes.
  # Example: assets/homo_sapiens_autosomes.txt

--archr_scrnaseq                    [false|path_to_matching_RNAseq_Seurat_object] Whether or not to integrate scRNAseq data.
--archr_scrnaseq_grouplist          [''|'example see conf/test.config']  scRNAseq cluster grouping info for constrained integration.
--custom_peaks                      [false|'example see conf/test.config']  For motif enrichment/deviation.

# Use below to exclude clusters for downstream analysis, refer to archr_clustering/Cluster_xxx_matrix.csv for valid cluster names
--filter_seurat_ilsi                [false|String]  Clusters to exclude for downstream analysis.
--filter_seurat_harmony             [false|String]  Clusters to exclude for downstream analysis.
  # Example: --filter_seurat_harmony 'C1, C2'

custom_peaks          = false // for motif enrichment/deviation module
  # Example: --custom_peaks 'Encode_K562_GATA1 = "https://www.encodeproject.org/files/ENCFF632NQI/@@download/ENCFF632NQI.bed.gz", Encode_GM12878_CEBPB = "https://www.encodeproject.org/files/ENCFF761MGJ/@@download/ENCFF761MGJ.bed.gz", Encode_K562_Ebf1 = "https://www.encodeproject.org/files/ENCFF868VSY/@@download/ENCFF868VSY.bed.gz", Encode_K562_Pax5 = "https://www.encodeproject.org/files/ENCFF339KUO/@@download/ENCFF339KUO.bed.gz"'
```

Module specific parameters can be adjusted by editing `conf/module.config`. Below are some examples, refer to `conf/module.config` for more examples.

 - For doublet detection with AMULET:
  ```
  'amulet_detect_doublets' {
    args = '--expectedoverlap 2 --maxinsertsize 2000'
  }
  ```
 - For plotting marker genes:
  ```
  'archr_marker_gene_clusters' {
    // args is for ArchR::getMarkerFeatures()
    args = 'useMatrix = "GeneScoreMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon"'

    // getMarkers_cutoff is for ArchR::getMarkers()
    getMarkers_cutoff = 'cutOff = "FDR <= 0.05 & Log2FC >= 1"'

    // marker_genes is for plotting marker genes
    marker_genes = 'default' // by default, the first 3 will be plotted, customized example see below:
    //marker_genes = 'CD34, GATA1, PAX5, MS4A1, EBF1, MME, CD14, CEBPB, MPO, IRF8, CD3D, CD8A, TBX21, IL7R'

    // args2 is for visualizing embedding
    args2 = 'colorBy = "GeneScoreMatrix", quantCut = c(0.01, 0.95)'

    // args3 is for track plotting with ArchR::ArchRBrowser()
    args3 = 'upstream = 50000, downstream = 50000'
  }
  ```

As shown above, for visualization purpose, scATACseq attempts to plot the first 3 marker genes (`marker_genes = 'default'`) from the resulting marker gene list. Such default parameters may not make direct sense to your research and are subject to user modifications. Below summarizes such parameters.
```
marker_genes from archr_marker_gene_clusters/clusters2 module: default to plot the first 3.
marker_genes from archr_scrnaseq_constrained module: default to plot the first 3.
marker_genes from archr_scrnaseq_unconstrained module: default to plot the first 3.
marker_genes from archr_marker_peaks_in_tracks_clusters/clusters2: default to plot the first 3.
marker_genes from archr_coaccessibility_clusters/clusters2: default to plot the first 3.
marker_genes from archr_peak2genelinkage_clusters2: default to plot the first 3.

use_groups/bgd_groups from archr_pairwise_test_clusters/clusters2: default to test the first group against the second.

motifs from archr_motif_deviations_clusters/clusters2: default to plot the first 3.
motifs from archr_footprinting_clusters/clusters2: default to plot the first 3.

trajectory_groups from archr_trajectory_clusters2: default to use the first 3 groups (since the order matters, this trajectory plotting is just for demo).
```

Note that modules related to clustering with only scATACseq data will be named as `xxx_clusters` whereas clustering based on integrated scRNAseq data will be named as `xxx_clusters2`.

## FASTQ files as input
FASTQ file paths (full path) must be saved into a **.csv** file (see below) and supplied with `--input_fastq`. Multiple runs of the same sample can be placed into separate rows with the same `sample_name`.

```
sample_name,path_fastq_1,path_fastq_2,path_barcode
SAMPLE_1,/full_path/xxx.fastq.gz,/full_path/xxx.fastq.gz,/full_path/xxx.fastq.gz
SAMPLE_2,/full_path/xxx.fastq.gz,/full_path/xxx.fastq.gz,/full_path/xxx.fastq.gz
```
| Column       | Description                                                                            |
|--------------|----------------------------------------------------------------------------------------|
| sample_name  | Must be identical for multiple runs of the same sample.                                |
| path_fastq_1 | Must use full path. File has to be gzipped and the extension of '.fq.gz' or 'fastq.gz'.|
| path_fastq_2 | Must use full path. File has to be gzipped and the extension of '.fq.gz' or 'fastq.gz'.|

An example .csv can be found [here](https://raw.githubusercontent.com/hukai916/scATACpipe/main/assets/example_samplesheet_fastq.csv).

In addition to the FASTQ files, you must also specify a preprocessing strategy with:
```
--preprocess    [default|10xgenomics|chromap] Preprocess strategy. Default to default.
```

The genome/annotation files are also required, and 3 options are available.

### Option 1: using UCSC/ENSEMBL genome
```
--ref_fasta_ensembl|--ref_fasta_ucsc    [string]  A genome name, either from ENSEMBL (e.g. homo_sapiens) or UCSC (e.g. mm10).
--species_latin_name                    [string]  Must be quoted, required if genome name not in ("hg38", "hg19", "mm10", "mm9").
```

### Option2: using custom genome
```
--ref_fasta             [string]  Path to refernce genome file.
--ref_gtf               [string]  Path to refernce gtf file.
--species_latin_name    [string]  Must be quoted.
```

### Option3: using existing genome index
If genome index files are readily available, you can skip the index-building step by directly supply the index folder.
```
--ref_bwa_index           [string]  Path to the bwa index folder. For '--preprocess default' only.
--ref_cellranger_index    [string]  Path to cellranger index folder. For '--preprocess 10xgenomics' only.
--ref_chromap_index       [string]  Path to chromap index folder. For '--preprocess chromap' only.
--species_latin_name    [string]  Must be quoted.
```

### Other parameters

Parameters related to downstream analysis can be tuned the same way as **Fragement files as input - Other parameters** above.

Similarly, **main pipeline parameters** must be supplied with command flags or configured inside `nextflow.confg`. These parameters are listed below (mainly relevant to PREPROCESS_DEFAULT).
```
--split_fastq           [true|false]  Whether or not to split fastq into chunks for more parallel jobs. Default to false.
--barcode_correction    [pheniqs|naive|false]  Barcode correction method, set to false to skip. Default to pheniqs.
--whitelist_barcode     [string]  Path to whitelist folder. Default to 'assets/whitelist_barcodes'.
--read1_adapter         [string]  For adapter trimming. Default to 'AGATCGGAAGAGC', which is Illumina standard adapters.
--read2_adapter         [string]  For adapter trimming. Default to 'AGATCGGAAGAGC', which is Illumina standard adapters.
--mapper                [bwa]  For mapping. Currently only support 'bwa'.
--filter                [both|improper|false]  For BAM file filtering: 'improper' means to filter out reads with low mapping quality or with extreme fragment size (outside of range 38 - 2000nt) etc. 'both' means to filter out both 'improper' and mitochondrial reads. Use false to skip filtering. Default to 'both'.
```

Also similarly, **module specific parameters** can be adjusted by editing corresponding sections in `conf/modules.config`. Below are some examples.
 - For barcode correction with pheniqs (PREPROCESS_DEFAULT):
  ```
  'correct_barcode_pheniqs' {
    read_count_cutoff = '10' // number of minimum reads to count as valid barcode
  }
  ```
 - For BAM file preparation (PREPROCESS_DEFAULT):
  ```
  prep_bam {
    args = '--shift_forward 4 --shift_reverse -5 --barcode_regex "[^:]*"'
  }
  ```
 - For tuning cellranger_atac_count (PREPROCESS_10XGENOMICS):
 ```
 'cellranger_atac_count' {
   args = '' // can be any natively support cellranger_atac_count flag
 }
 ```

## Running the pipeline
For example commands, see [Quick Start](https://github.com/hukai916/scATACpipe/#quick-start).

## Version control
To use a specific version of the pipeline:
```bash
cd scATACpipe
git pull # get all updates
git branch # show all available braches
git checkout dev # switch to dev branch
```

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/scatacpipe --input samplesheet.csv --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/scatacpipe
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/scatacpipe releases page](https://github.com/nf-core/scatacpipe/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `conda`
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: STAR_ALIGN {
        memory = 100.GB
    }
}
```

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Tool-specific options

For the ultimate flexibility, we have implemented and are using Nextflow DSL2 modules in a way where it is possible for both developers and users to change tool-specific command-line arguments (e.g. providing an additional command-line argument to the `STAR_ALIGN` process) as well as publishing options (e.g. saving files produced by the `STAR_ALIGN` process that aren't saved by default by the pipeline). In the majority of instances, as a user you won't have to change the default options set by the pipeline developer(s), however, there may be edge cases where creating a simple custom config file can improve the behaviour of the pipeline if for example it is failing due to a weird error that requires setting a tool-specific parameter to deal with smaller / larger genomes.

The command-line arguments passed to STAR in the `STAR_ALIGN` module are a combination of:

* Mandatory arguments or those that need to be evaluated within the scope of the module, as supplied in the [`script`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L49-L55) section of the module file.

* An [`options.args`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L56) string of non-mandatory parameters that is set to be empty by default in the module but can be overwritten when including the module in the sub-workflow / workflow context via the `addParams` Nextflow option.

The nf-core/rnaseq pipeline has a sub-workflow (see [terminology](https://github.com/nf-core/modules#terminology)) specifically to align reads with STAR and to sort, index and generate some basic stats on the resulting BAM files using SAMtools. At the top of this file we import the `STAR_ALIGN` module via the Nextflow [`include`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/subworkflows/nf-core/align_star.nf#L10) keyword and by default the options passed to the module via the `addParams` option are set as an empty Groovy map [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/subworkflows/nf-core/align_star.nf#L5); this in turn means `options.args` will be set to empty by default in the module file too. This is an intentional design choice and allows us to implement well-written sub-workflows composed of a chain of tools that by default run with the bare minimum parameter set for any given tool in order to make it much easier to share across pipelines and to provide the flexibility for users and developers to customise any non-mandatory arguments.

When including the sub-workflow above in the main pipeline workflow we use the same `include` statement, however, we now have the ability to overwrite options for each of the tools in the sub-workflow including the [`align_options`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/workflows/rnaseq.nf#L225) variable that will be used specifically to overwrite the optional arguments passed to the `STAR_ALIGN` module. In this case, the options to be provided to `STAR_ALIGN` have been assigned sensible defaults by the developer(s) in the pipeline's [`modules.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L70-L74) and can be accessed and customised in the [workflow context](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/workflows/rnaseq.nf#L201-L204) too before eventually passing them to the sub-workflow as a Groovy map called `star_align_options`. These options will then be propagated from `workflow -> sub-workflow -> module`.

As mentioned at the beginning of this section it may also be necessary for users to overwrite the options passed to modules to be able to customise specific aspects of the way in which a particular tool is executed by the pipeline. Given that all of the default module options are stored in the pipeline's `modules.config` as a [`params` variable](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L24-L25) it is also possible to overwrite any of these options via a custom config file.

Say for example we want to append an additional, non-mandatory parameter (i.e. `--outFilterMismatchNmax 16`) to the arguments passed to the `STAR_ALIGN` module. Firstly, we need to copy across the default `args` specified in the [`modules.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L71) and create a custom config file that is a composite of the default `args` as well as the additional options you would like to provide. This is very important because Nextflow will overwrite the default value of `args` that you provide via the custom config.

As you will see in the example below, we have:

* appended `--outFilterMismatchNmax 16` to the default `args` used by the module.
* changed the default `publish_dir` value to where the files will eventually be published in the main results directory.
* appended `'bam':''` to the default value of `publish_files` so that the BAM files generated by the process will also be saved in the top-level results directory for the module. Note: `'out':'log'` means any file/directory ending in `out` will now be saved in a separate directory called `my_star_directory/log/`.

```nextflow
params {
    modules {
        'star_align' {
            args          = "--quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --readFilesCommand zcat --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outFilterMismatchNmax 16"
            publish_dir   = "my_star_directory"
            publish_files = ['out':'log', 'tab':'log', 'bam':'']
        }
    }
}
```

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

    * For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Singularity:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
