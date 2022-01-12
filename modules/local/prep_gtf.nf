// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREP_GTF {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir: 'prep_gtf', publish_id:'') }
    container "hukai916/miniconda3_bio:0.1"

    input:
    path genome_fasta
    val genome_name
    path gtf

    output:
    path "final.gtf.gz", emit: gtf
    val genome_name, emit: genome_name

    script:

    """
    if [[ $gtf == *.gz ]]; then
      gunzip -c $gtf > annotation.gtf
    else
      if [[ $gtf != annotation.gtf ]]; then
        mv $gtf annotation.gtf
      fi
    fi

    # convert to GTF in case GFF3 is supplied, which seems not compatible with ArchR
    gffread annotation.gtf -T -o- > annotation.v2.gtf

    # sort the gtf by gene_id in case some entries that belong to the same gene are not in order
    sort_gtf.py annotation.v2.gtf > sorted.gtf

    # add 'chr' prefix: required by ArchR
    cat sorted.gtf | awk 'BEGIN{FS=OFS="\\t"}/^#/{print; next} NF==9{\$1 = gensub(/(chr)?(.+).*/, "chr\\\\2", "g", \$1); print}' > chrPrefixed.sorted.gtf

    # make sure gtf config is a subset of genome.fa config, required for cellranger_index step
    extract_gtf.py $genome_fasta chrPrefixed.sorted.gtf > subset.chrPrefixed.sorted.gtf

    # make sure feature end coordinates don't exceed config boundaries, required for cellranger_index step
    filter_gtf.py $genome_fasta subset.chrPrefixed.sorted.gtf > filtered.subset.chrPrefixed.sorted.gtf

    # add 'gene' entries in case not there
    ## step1: add 'gene' with gffread, output to gff3
    gffread -E --keep-genes filtered.subset.chrPrefixed.sorted.gtf -o- > final.gff3
    ## step2: convert back to gtf with bioinfokit since the conversion using gffread is buggy
    gff3_to_gtf.py final.gff3

    gzip final.gtf
    """
}
