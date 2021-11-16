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
        def command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
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
            log.info "${colors.cyan}${"\nSupplied genome not supported! Try '--support_genome'.\n"}${colors.reset}\n\n"
            System.exit(0)

          }
        } else if (params.ref_fasta_ucsc) {
          def res = NfcoreSchema.checkGenome(workflow, params).toBoolean()
          if (!res) {
            Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
            log.info "${colors.cyan}${"\nSupplied genome not supported! Try '--support_genome'.\n"}${colors.reset}\n\n"
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

        // Check --input_archr/--input_preprocess has been provided
        if ((!params.input_preprocess) && (!params.input_archr)) {
            println ''
            log.error "Pls supply either fragment file (--input_archr) or samplesheet tsv file (--input_preprocess)!"
            println ''
            System.exit(0)
        } else if (params.input_archr) {
          // 3 DOWNSTREAM options
            if ((!params.archr_genome) && (!params.archr_genome_fasta || !params.archr_gtf || !params.species_latin_name) && (!params.archr_bsgenome || !params.archr_txdb || !params.archr_org)) {
              println ''
              def out_string = "Insufficient parameters supplied for DOWNSTREAM analysis!\nPls choose one:\n\n"
              out_string += "Option1:\n  --archr_genome [a genome name]\n  --archr_blacklist [optional, path to blacklist file]\n"
              out_string += "Option2:\n  --archr_genome_fasta [path to genome fasta]\n  --archr_gtf [path to gtf file]\n  --species_latin_name [latin name for genome organism, must be quoted]\n  --archr_blacklist [optional, path to blacklist file]\n"
              out_string += "Option3:\n  --archr_bsgenome [Bioconductor BSgenome name]\n  --archr_txdb [Bioconductor TxDb name]\n  --archr_org [Bioconductor OrgDb name]\n  --archr_blacklist [optional, path to blacklist file]"

              log.error out_string
              println ''
              System.exit(0)
            } else {
            log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
            }
        } else if (params.input_preprocess) {
          // First check if PREPROCESS parameters satisfied
          if (params.preprocess == "default") {
            if ((!params.ref_fasta_ucsc) && (!params.ref_fasta_ensembl) && (!params.ref_fasta || !params.ref_gtf) && (!params.ref_bwa_index || !params.ref_minimap2_index)) {
              println ''
              def out_string = "Insufficient parameters supplied for PREPROCESS_DEFAULT!\n\n"
              out_string += "Option1 (index folder):\n  --ref_bwa_index | --ref_minimap2_index [path to index folder]\n"
              out_string += "Option2 (custom genome fasta & gtf):\n  --ref_fasta [path to genome fasta file]\n  --ref_gtf [path to gtf file]\n"
              out_string += "Option3 (UCSC genome):\n  --ref_fasta_ucsc [UCSC genome]\n"
              out_string += "Option4 (ENSEMBL genome):\n  --ref_fasta_ensembl [ENSEMBL genome name]\n\n"
              log.error out_string
              System.exit(0)
            } else {
              log.info "Validating input params for DOWNSTREAM_DEFAULT, passed."
            }
          } else if (params.preprocess == "10xgenomics") {
            if ((!params.ref_fasta_ucsc) && (!params.ref_fasta_ensembl) && (!params.ref_fasta || !params.ref_gtf) && (!params.ref_cellranger_index)) {
              println ''
              def out_string = "Insufficient parameters supplied for PREPROCESS_10XGENOMICS!\n\n"
              out_string += "Option1 (index folder):\n  --ref_cellranger_index [path to index folder]\n"
              out_string += "Option2 (custom genome fasta & gtf):\n  --ref_fasta [path to genome fasta file]\n  --ref_gtf [path to gtf file]\n"
              out_string += "Option3 (UCSC genome):\n  --ref_fasta_ucsc [UCSC genome]\n"
              out_string += "Option4 (ENSEMBL genome):\n  --ref_fasta_ensembl [ENSEMBL genome name]\n\n"
              log.error out_string
              System.exit(0)
            } else {
              log.info "Validating input params for DOWNSTREAM_10XGENOMICS, passed."
            }
          } else {
            log.error "Pls supply --preprocess [default | 10xgenomics]"
            System.exit(0)
          }

          // Also check if DOWNSTREAM_ARCHR parameters satisfied
          if (params.preprocess == "default") {
            if (params.ref_bwa_index || params.ref_minimap2_index) {
              if (!params.ref_fasta_ensembl && !params.ref_fasta_ucsc) {
                log.error "Pls supply --ref_fasta_ensembl [ENSEMBL genome name] | --ref_fasta_ucsc [UCSC genome name]"
                System.exit(0)
              } else {
                log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
              }
            } else if (params.ref_fasta) {
              if (!params.archr_gtf) {
                log.error "Pls supply --archr_gtf [path to gtf file]"
                System.exit(0)
              } else {
                log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
              }
            } else if (params.ref_fasta_ensembl || params.ref_fasta_ucsc) {
              log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
            }
          } else if (params.preprocess == "10xgenomics") {
              if (params.ref_cellranger_index) {
                if (!params.ref_fasta_ensembl && !params.ref_fasta_ucsc) {
                  log.error "Pls supply --ref_fasta_ensembl [ENSEMBL genome name] | --ref_fasta_ucsc [UCSC genome name]"
                  System.exit(0)
                } else {
                  log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
                }
              } else if (params.ref_fasta) {
                if (!params.archr_gtf) {
                  log.error "Pls supply --archr_gtf [path to gtf file]"
                  System.exit(0)
                } else {
                  log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
                }
              } else if (params.ref_fasta_ensembl || params.ref_fasta_ucsc) {
                log.info "Validating input params for DOWNSTREAM_ARCHR, passed."
              }
          }
        }
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
