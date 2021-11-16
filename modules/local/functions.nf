/*
 * -----------------------------------------------------
 *  Utility functions used in nf-core DSL2 module files
 * -----------------------------------------------------
 */

/*
 * Extract name of software tool from process name using $task.process
 */
def getSoftwareName(task_process) {
    return task_process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()
}

/*
 * Function to initialise default values and to generate a Groovy Map of available options for nf-core modules
 */
def initOptions(Map args) {
    def Map options = [:]
    options.args          = args.args ?: ''
    options.args2         = args.args2 ?: ''
    options.args3         = args.args3 ?: ''
    options.read_count_cutoff = args.read_count_cutoff ?: ''
    options.marker_genes  = args.marker_genes ?: ''
    options.cutoff        = args.cutoff ?: ''
    options.motifs        = args.motifs ?: ''
    options.norm_method   = args.norm_method ?: ''
    options.tss_norm_method = args.tss_norm_method ?: ''
    options.tss_flank     = args.tss_flank ?: ''
    options.flank_norm    = args.flank_norm ?: ''
    options.bam_coverage  = args.bam_coverage ?: ''
    options.gene_to_color = args.gene_to_color ?: ''
    options.publish_by_id = args.publish_by_id ?: false
    options.publish_dir   = args.publish_dir ?: ''
    options.publish_files = args.publish_files
    options.suffix        = args.suffix ?: ''
    return options
}

/*
 * Tidy up and join elements of a list to return a path string
 */
def getPathFromList(path_list) {
    def paths = path_list.findAll { item -> !item?.trim().isEmpty() }  // Remove empty entries
    paths = paths.collect { it.trim().replaceAll("^[/]+|[/]+\$", "") } // Trim whitespace and trailing slashes
    return paths.join('/')
}

/*
 * Function to save/publish module results
 */
def saveFiles(Map args) {
    if (!args.filename.endsWith('.version.txt')) {
        def ioptions = initOptions(args.options)
        def path_list = [ ioptions.publish_dir ?: args.publish_dir ]
        if (ioptions.publish_by_id) {
            path_list.add(args.publish_id)
        }
        if (ioptions.publish_files instanceof Map) {
            for (ext in ioptions.publish_files) {
                if (args.filename.endsWith(ext.key)) {
                    def ext_list = path_list.collect()
                    ext_list.add(ext.value)
                    return "${getPathFromList(ext_list)}/$args.filename"
                }
            }
        } else if (ioptions.publish_files == null) {
            return "${getPathFromList(path_list)}/$args.filename"
        }
    }
}

// Function to prepare for ArchR genome name:
// def get_bsgenome(archr_genome, archr_custom_genome, archr_txdb, archr_org, archr_bsgenome, ref_fasta_ucsc, ref_fasta_ensembl, ref_cellranger_ucsc, ref_cellranger_ensembl) {
//   // Natively supported ArchR genomes:
//   def archr_support_genome = ["hg19", "hg38", "mm9", "mm10"]
//   // Other supported ArchR bsgenomes with ArchR create annotation functions:
//   def archr_custom_bsgenome = ["bosTau9", "ce11", "canFam3", "danRer11", "dm6", "galGal6", "rheMac10", "panTro6", "rn6", "sacCer3", "susScr11"]
//
//   if (archr_custom_genome == "yes") {
//     // Check if txdb, org, and bsgenome are specified
//     if (!(archr_txdb) || !(archr_org) || !(archr_bsgenome)) {
//       exit 1, '--archr_custom_genome set to "yes", pls also supply --archr_txdb, --archr_org, and --archr_bsgenome.'
//     } else {
//       return ["custom", "custom"]
//     }
//   } else if (archr_custom_genome == "no") {
//     if (archr_support_genome.contains(params.archr_genome)) {
//       return [archr_genome, "ready"]
//     }
//     if (archr_custom_bsgenome.contains(params.archr_genome)) {
//       return [archr_genome, "need_build"]
//     }
//     if (archr_support_genome.contains(ref_fasta_ucsc)) {
//       return [ref_fasta_ucsc, "ready_ucsc"]
//     }
//     if (archr_support_genome.contains(ref_cellranger_ucsc)) {
//       return [ref_cellranger_ucsc, "ready_ucsc"]
//     }
//     if (archr_custom_bsgenome.contains(ref_fasta_ucsc)) {
//       return [ref_fasta_ucsc, "need_build_ucsc"]
//     }
//     if (archr_support_genome.contains(ensembl2ucsc[ref_fasta_ensembl])) {
//       return [ensembl2ucsc[ref_fasta_ensembl], "ready_ensembl"]
//     }
//     if (archr_support_genome.contains(ensembl2ucsc[ref_cellranger_ensembl])) {
//       return [ensembl2ucsc[ref_cellranger_ensembl], "ready_ensembl"]
//     }
//     if (archr_custom_bsgenome.contains(ref_fasta_ensembl)) {
//       return [ref_fasta_ucsc, "need_build_ensembl"]
//     }
//   } else {
//     exit 1, '--archr_custom_genome must be either "yes" or "no".'
//   }
//
//   return ["NA", "NA"]
// }



// // Map: ensemble name to filename:
// def get_ensembl_filename() {
//       def dict_genome_name = [homo_sapiens: "Homo_sapiens.GRCh38.dna.primary_assembly", mus_musculus: "Mus_musculus.GRCm39.dna.primary_assembly", bos_taurus: "Bos_taurus.ARS-UCD1.2.dna.toplevel", caenorhabditis_elegans: "Caenorhabditis_elegans.WBcel235.dna.toplevel", danio_rerio: "Danio_rerio.GRCz11.dna.primary_assembly", drosophila_melanogaster: "Drosophila_melanogaster.BDGP6.32.dna.toplevel", gallus_gallus: "Gallus_gallus.GRCg6a.dna.toplevel", macaca_mulatta: "Macaca_mulatta.Mmul_10.dna.toplevel",  pan_troglodytes: "Pan_troglodytes.Pan_tro_3.0.dna.toplevel", rattus_norvegicus: "Rattus_norvegicus.Rnor_6.0.dna.toplevel", saccharomyces_cerevisiae: "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel",  sus_scrofa: "Sus_scrofa.Sscrofa11.1.dna.toplevel"
//       ]
//       return dict_genome_name
// }
