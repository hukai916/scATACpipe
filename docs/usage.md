# scATACpipe: Usage
(Some parts adapted from nf-core [TEMPLATE](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/docs/usage.md).)

## Table of Contents
[Introduction](#introduction)

[Basics](#basics)

[Fragment files as input](#fragment-files-as-input)
  * [Option1: using UCSC/ENSEMBL genome](#option1-using-ucscensembl-genome)
  * [Option2: using custom genome](#option2-using-custom-genome)
  * [Option3: using Bioconductor annotations](#option3-using-bioconductor-annotations)
  * [Other parameters](#other-parameters)

[FASTQ files as input](#fastq-files-as-input)
  * [Option1: using UCSC/ENSEMBL genome](#option1-using-ucscensembl-genome-1)
  * [Option2: using custom genome](#option2-using-custom-genome-1)
  * [Option3: using existing genome index](#option3-using-existing-genome-index)
  * [Other parameters](#other-parameters-1)

[Running the pipeline](#running-the-pipeline)

[Core Nextflow arguments](#core-nextflow-arguments)
  * [-profile](#-profile)
  * [-resume](#-resume)
  * [-c](#-c)

[Custom configuration](#custom-configuration)
  * [Resource requests](#resource-requests)
  * [Module-specific options](#module-specific-options)
  * [Create your own config files](#create-your-own-config-files)

[Running in the background](#running-in-the-background)

[Nextflow memory requirements](#nextflow-memory-requirements)

[Version control](#version-control)

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

### Option1: using UCSC/ENSEMBL genome
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
  ```nextflow
  'amulet_detect_doublets' {
    args = '--expectedoverlap 2 --maxinsertsize 2000'
  }
  ```
 - For plotting marker genes:
  ```nextflow
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

### Option1: using UCSC/ENSEMBL genome
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
--species_latin_name      [string]  Must be quoted.
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
  ```nextflow
  'correct_barcode_pheniqs' {
    read_count_cutoff = '10' // number of minimum reads to count as valid barcode
  }
  ```
 - For BAM file preparation (PREPROCESS_DEFAULT):
  ```nextflow
  prep_bam {
    args = '--shift_forward 4 --shift_reverse -5 --barcode_regex "[^:]*"'
  }
  ```
 - For tuning cellranger_atac_count (PREPROCESS_10XGENOMICS):
  ```nextflow
  'cellranger_atac_count' {
   args = '' // can be any natively support cellranger_atac_count flag
  }
  ```

## Running the pipeline
For example commands, see [Quick Start](https://github.com/hukai916/scATACpipe/#quick-start).

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### -profile
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

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

#### -resume
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### -c
Specify the path to a specific config file. See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests
Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/hukai916/scATACpipe/blob/dev/conf/base.config/#L21) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the scATACpipe is failing after multiple re-submissions of the `FASTQC` process due to an exit code of `137` this would indicate that there is an out of memory issue.

To bypass this error you would need to find exactly which resources are set by the `FASTQC` process. The quickest way is to search for `process FASTQC` in the [hukai916/scATACpipe GitHub repo](https://github.com/hukai916/scATACpipe/search?q=process+FASTQC). We have standardized the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/local/fastqc.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_low`](https://github.com/hukai916/scATACpipe/blob/93762729a65f96d831bed9ccf97406ed2e3b4166/modules/local/fastqc.nf#L8). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organize workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_low` label are set in the pipeline's [`base.config`](https://github.com/hukai916/scATACpipe/blob/93762729a65f96d831bed9ccf97406ed2e3b4166/conf/base.config#L38-L42) which in this case is defined as 12GB. Providing you haven't set any other standard Nextflow parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `FASTQC` process failure by creating a custom config file that sets at least 12GB of memory, in this case increased to 20GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: FASTQC {
        memory = 20.GB
    }
}
```

> **NB:** We specify just the process name i.e. `FASTQC` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running.

### Module-specific options
For the ultimate flexibility, we have implemented and are using Nextflow DSL2 modules in a way where it is possible for both developers and users to change module-specific (each module usually wraps around one or more tools) command-line arguments (e.g. providing an additional command-line argument to the `FASTQC` process) as well as publishing options (e.g. saving files produced by certain process that aren't saved by default by the pipeline).

The command-line arguments passed to FASTQC in the `FASTQC` module are a combination of:

* Mandatory arguments or those that need to be evaluated within the scope of the module, as supplied in the [`fastqc.nf`](https://github.com/hukai916/scATACpipe/blob/main/modules/local/fastqc.nf#L27-L30) section of the module file.

* An [`options.args`](https://github.com/hukai916/scATACpipe/blob/main/modules/local/fastqc.nf#L26) string of non-mandatory parameters that is set to be empty by default in the module but can be overwritten when including the module in the sub-workflow / workflow context via the `addParams` Nextflow option.

Basically, the sub-workflow __PREPROCESS_DEFAULT__ imports the `FASTQC` module via the Nextflow [`include`](https://github.com/hukai916/scATACpipe/blob/main/workflows/preprocess_default.nf#L42) keyword and by default the options passed to the module via the `addParams` option are set as an empty Groovy map [here](https://github.com/hukai916/scATACpipe/blob/main/modules/local/functions.nf#L18-L19); this in turn means `options.args` will be set to empty by default in the module file too. This is an intentional design choice and allows us to implement well-written sub-workflows composed of a chain of tools that by default run with the bare minimum parameter set for any given tool in order to make it much easier to share across pipelines and to provide the flexibility for users and developers to customize any non-mandatory arguments.

When including the sub-workflow above in the main pipeline workflow we use the same `include` statement, however, we now have the ability to overwrite options for each of the tools in the sub-workflow. In this case, the options to be provided to `FASTQC` have been assigned sensible defaults by the developer(s) in the pipeline's [`conf/modules.config`](https://github.com/hukai916/scATACpipe/blob/main/conf/modules.config#L21-L23), which can be accessed and customized there. Alternatively, these parameters can be overwritten with a custom config file as described in [-c](#-c).

Say for example we want to append an additional, non-mandatory parameter (i.e. `--noextract`) to the arguments passed to the `FASTQC` module. Firstly, we need to copy across the default `args` specified in the [`conf/modules.config`](https://github.com/hukai916/scATACpipe/blob/main/conf/modules.config#L22) and create a custom config file that is a composite of the default `args` as well as the additional options you would like to provide. This is very important because Nextflow will overwrite the default value of `args` that you provide via the custom config.

As you will see in the example below, we have:

* appended `noextract` to the default `args (--quiet)` used by the module.
* changed the default `publish_dir` value to where the files will eventually be published in the main results directory.

```nextflow
params {
    modules {
        'fastqc' {
          args = '--quiet --noextract'
        }
    }
}
```

### Create your own config files
See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

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

## Version control
By default, pipeline from the `main` branch will be used, which is identical to the latest released version.

To use a specific version of the pipeline:
```bash
cd scATACpipe
git pull # get all updates
git branch -r # show all available braches
git checkout origin/dev # switch to dev branch
```
