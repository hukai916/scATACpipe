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

    # convert to GTF in case GFF3 is supplied: ArchR seems only allow GTF
    gffread annotation.gtf -T -o- > annotation.v2.gtf

    # add 'chr' prefix: required by ArchR
    cat annotation.v2.gtf | awk 'BEGIN{FS=OFS="\\t"}/^#/{print; next} NF==9{\$1 = gensub(/(chr)?(.+).*/, "chr\\\\2", "g", \$1); print}' > chrPrefixed.gtf

    # make sure gtf config is a subset of genome.fa config, required for cellranger_index step
    extract_gtf.py $genome_fasta chrPrefixed.gtf > subset.chrPrefixed.gtf

    # sort the gtf by gene_id in case some entries that belong to the same gene are not in order, also output gene_ranges.tsv
    sort_gtf.py subset.chrPrefixed.gtf > sorted.subset.chrPrefixed.gtf

    # add 'gene' entries in case not there
      ## step1: add 'gene' with gffread, output to gff3
    gffread -E --keep-genes sorted.subset.chrPrefixed.gtf -o- > sorted.subset.chrPrefixed.gff3
      ## step2: correct the 'gene range' added by gffread: psudo-autosomals: https://github.com/gpertea/gffread/issues/86
    correct_gtf_gene_range.py sorted.subset.chrPrefixed.gff3 gene_ranges.tsv > final.gff3
      ## step3: convert back to gtf with bioinfokit since the conversion using gffread is buggy
    gff3_to_gtf.py final.gff3

    gzip final.gtf

    # clean-up
    rm annotation.gtf annotation.v2.gtf chrPrefixed.sorted.gtf final.gff3 sorted.gtf subset.chrPrefixed.sorted.gff3 subset.chrPrefixed.sorted.gtf gene_ranges.tsv
    """
}
