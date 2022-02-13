//
// This file holds several functions specific to the main.nf workflow in the nf-core/scatacpipe pipeline
//

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {

        def command1 = "nextflow run main.nf --input_fragment samplesheet.csv --archr_genome hg19 -profile [docker|singularity]"
        def command2 = "nextflow run main.nf --input_fastq samplesheet.csv --preprocess [default|10xgenomics] --ref_fasta_ucsc hg19 --species_latin_name 'homo sapiens' -profile [docker|singularity]"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command1, command2)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Print supported genome to screen if required
    //
    public static String support_genome(workflow, params, log) {
        def command = "nextflow run main.nf --support_genome"
        def genome_string = ''
        genome_string += NfcoreSchema.paramsGenome(workflow, params, command)
        genome_string += '\n'
        genome_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return genome_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Print supported genome to screen if required
        if (params.support_genome) {
          log.info support_genome(workflow, params, log)
          System.exit(0)
        }

        // Prompt if supplied genome not supported
        if (params.ref_fasta_ensembl) {
          def res = NfcoreSchema.checkGenome(workflow, params).toBoolean()
          if (!res) {
            Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
            log.error "Supplied genome not supported! Try '--support_genome'.\n"
            System.exit(0)
          }
        } else if (params.ref_fasta_ucsc) {
          def res = NfcoreSchema.checkGenome(workflow, params).toBoolean()
          if (!res) {
            Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
            log.error "Supplied genome not supported! Try '--support_genome'.\n"
            System.exit(0)
          }
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check the hostnames against configured profiles
        NfcoreTemplate.hostName(workflow, params, log)

        // Check --input_fragment/--input_fastq has been provided
        if ((!params.input_fastq) && (!params.input_fragment)) {
            println ''
            log.error "Pls supply either fragment (--input_fragment) or fastq samplesheet file (--input_fastq)!"
            println ''
            System.exit(0)
        } else if (params.input_fragment) {
          // 3 DOWNSTREAM options
            // Check if DOWNSTREAM parameters satisfied
            log.info "Validating input params for DOWNSTREAM_ARCHR ..."
            if ((!params.archr_genome) && (!params.archr_genome_fasta || !params.ref_gtf || !params.species_latin_name) && (!params.archr_bsgenome || !params.archr_txdb || !params.archr_org)) {
              println ''
              def out_string = "Insufficient parameters supplied for DOWNSTREAM analysis!\nPls choose one:\n\n"
              out_string += "Option1:\n  --archr_genome [a genome name]\n  --archr_blacklist [optional, path to blacklist file]\n"
              out_string += "Option2:\n  --archr_genome_fasta [path to genome fasta]\n  --ref_gtf [path to gtf file]\n  --species_latin_name [latin name for genome organism, must be quoted]\n  --archr_blacklist [optional, path to blacklist file]\n"
              out_string += "Option3:\n  --archr_bsgenome [Bioconductor BSgenome name]\n  --archr_txdb [Bioconductor TxDb name]\n  --archr_org [Bioconductor OrgDb name]\n  --archr_blacklist [optional, path to blacklist file]"

              log.error out_string
              println ''
              System.exit(0)
            }
            log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
        } else if (params.input_fastq) {
          // Check if PREPROCESS parameters satisfied
          if (params.preprocess == "default") {
            log.info "Validating input params for PREPROCESS_DEFAULT ..."
            if ((!params.ref_fasta_ucsc) && (!params.ref_fasta_ensembl) && (!params.ref_fasta || !params.ref_gtf) && (!params.ref_bwa_index)) {
              println ''
              def out_string = "Insufficient parameters supplied for PREPROCESS_DEFAULT!\n\n"
              out_string += "Option1 (index folder):\n  --ref_bwa_index [path to index folder]\n"
              out_string += "Option2 (custom genome fasta & gtf):\n  --ref_fasta [path to genome fasta file]\n  --ref_gtf [path to gtf file]\n"
              out_string += "Option3 (UCSC genome):\n  --ref_fasta_ucsc [UCSC genome]\n"
              out_string += "Option4 (ENSEMBL genome):\n  --ref_fasta_ensembl [ENSEMBL genome name]\n\n"
              log.error out_string
              System.exit(0)
            }
            log.info "Validating input params for PREPROCESS_DEFAULT, passed."
          } else if (params.preprocess == "10xgenomics") {
            log.info "Validating input params for PREPROCESS_10XGENOMICS ..."
            if ((!params.ref_fasta_ucsc) && (!params.ref_fasta_ensembl) && (!params.ref_fasta || !params.ref_gtf) && (!params.ref_cellranger_index)) {
              println ''
              def out_string = "Insufficient parameters supplied for PREPROCESS_10XGENOMICS!\n\n"
              out_string += "Option1 (index folder):\n  --ref_cellranger_index [path to index folder]\n"
              out_string += "Option2 (custom genome fasta & gtf):\n  --ref_fasta [path to genome fasta file]\n  --ref_gtf [path to gtf file]\n"
              out_string += "Option3 (UCSC genome):\n  --ref_fasta_ucsc [UCSC genome]\n"
              out_string += "Option4 (ENSEMBL genome):\n  --ref_fasta_ensembl [ENSEMBL genome name]\n\n"
              log.error out_string
              System.exit(0)
            }
            log.info "Validating input params for PREPROCESS_10XGENOMICS, passed."
          } else if (params.preprocess == "chromap") {
            log.info "Validating input params for PREPROCESS_CHROMAP ..."
            if ((!params.ref_fasta_ucsc) && (!params.ref_fasta_ensembl) && (!params.ref_fasta || !params.ref_gtf)) {
              println ''
              def out_string = "Insufficient parameters supplied for PREPROCESS_CHROMAP!\n\n"
              out_string += "Option1 (custom genome fasta & gtf):\n  --ref_fasta [path to genome fasta file]\n  --ref_gtf [path to gtf file]\n"
              out_string += "Option2 (UCSC genome):\n  --ref_fasta_ucsc [UCSC genome]\n"
              out_string += "Option3 (ENSEMBL genome):\n  --ref_fasta_ensembl [ENSEMBL genome name]\n"
              out_string += "For all options, you can supply '--ref_chromap_index [path to chromap index]' to skip building index.\n\n"
              log.error out_string
              System.exit(0)
            }
            log.info "Validating input params for PREPROCESS_CHROMAP, passed."
          } else {
            log.error "Pls supply --preprocess [default | 10xgenomics | chromap]"
            System.exit(0)
          }

          // Also check if DOWNSTREAM_ARCHR parameters satisfied
          log.info "Validating input params for DOWNSTREAM_ARCHR ..."
          if (params.preprocess == "default") {
            if (params.ref_bwa_index) {
              if (!(params.ref_fasta_ensembl && params.species_latin_name) && !(params.ref_fasta_ucsc && params.species_latin_name)) {
                log.error "Pls supply --ref_fasta_ensembl [ENSEMBL genome name] | --ref_fasta_ucsc [UCSC genome name]\nPls also supply --species_latin_name [Must be quoted]"
                System.exit(0)
              }
            } else if (params.ref_fasta) {
              if (!params.ref_gtf || !params.species_latin_name) {
                log.error "Pls supply --ref_gtf [path to gtf file] AND --species_latin_name [Must be quoted]"
                System.exit(0)
              }
            } else if (params.ref_fasta_ensembl || params.ref_fasta_ucsc) {
              if (!params.species_latin_name) {
                log.error "Pls also supply --species_latin_name [Must be quoted]"
                System.exit(0)
              }
            }
          } else if (params.preprocess == "10xgenomics") {
              if (params.ref_cellranger_index) {
                if (!(params.ref_fasta_ensembl && params.species_latin_name) && !(params.ref_fasta_ucsc && params.species_latin_name)) {
                  log.error "Pls supply --ref_fasta_ensembl [ENSEMBL genome name] | --ref_fasta_ucsc [UCSC genome name]\nPls also supply --species_latin_name [Must be quoted]"
                  System.exit(0)
                }
              } else if (params.ref_fasta) {
                if (!params.ref_gtf || !params.species_latin_name) {
                  log.error "Pls supply --ref_gtf [path to gtf file] AND --species_latin_name [Must be quoted]"
                  System.exit(0)
                }
              } else if (params.ref_fasta_ensembl || params.ref_fasta_ucsc) {
                if (!params.species_latin_name) {
                  log.error "Pls also supply --species_latin_name [Must be quoted]"
                  System.exit(0)
                }
              }
          } else if (params.preprocess == "chromap") {
            if (params.ref_fasta) {
              if (!params.ref_gtf || !params.species_latin_name) {
                log.error "Pls supply --ref_gtf [path to gtf file] AND --species_latin_name [Must be quoted]"
                System.exit(0)
              }
            } else if (params.ref_fasta_ensembl || params.ref_fasta_ucsc) {
              if (!params.species_latin_name) {
                log.error "Pls also supply --species_latin_name [Must be quoted]"
                System.exit(0)
              }
            }
          }
          log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
        }

        // Check if other parameters are acceptable
        log.info "Validating other parameters ..."
        if (params.mapper) {
          if (!(params.mapper == "bwa") && !(params.mapper == "bowtie2")) {
            log.error "--mapper must be from 'bwa', 'bowtie2(todo)'."
            System.exit(0)
          }
        }
        if (params.barcode_correction) {
          if (!(params.barcode_correction == "naive") && !(params.barcode_correction == "pheniqs")) {
            log.error "--barcode_correction must be from 'naive', 'pheniqs'."
            System.exit(0)
          }
        }
        if (params.filter) {
          if (!(params.filter == "both") && !(params.filter == "improper")) {
            log.error "--filter must be from 'both', 'improper'."
            System.exit(0)
          }
        }
        if (params.doublet_removal_algorithm) {
          if (!(params.doublet_removal_algorithm == 'archr') && !(params.doublet_removal_algorithm == 'amulet')) {
            log.error "--doublet_removal_algorithm must be from 'archr', 'amulet', or false."
            System.exit(0)
          } else if (params.doublet_removal_algorithm == 'archr') {
            if (!params.archr_filter_doublets_ratio) {
              log.error "--archr_filter_doublets_ratio must be supplied!"
              System.exit(0)
            }
          } else if (params.doublet_removal_algorithm == 'amulet') {
            if (!(params.amulet_rmsk_bed) || !(params.amulet_autosomes)) {
              log.error "both --amulet_rmsk_bed and --amulet_autosomes must be supplied!"
              System.exit(0)
            }
          }
        }
        log.info "Validating other parameters, passed."
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttribute(params, attribute) {
        def val = ''
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
    }
}
