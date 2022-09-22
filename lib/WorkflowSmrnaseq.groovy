//
// This file holds several functions specific to the workflow/smrnaseq.nf in the nf-core/smrnaseq pipeline
//

class WorkflowSmrnaseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        // if (!params.fasta) {
        //     log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        //     System.exit(1)
        // }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }//
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

    /*
    * Format the protocol
    * Given the protocol parameter (params.protocol),
    * this function formats the protocol such that it is fit for the respective
    * subworkflow
    */
    public static void formatProtocol(protocol) {

    switch(protocol){
        case 'illumina':
            params.replace("clip_r1", 0)
            params.replace("three_prime_clip_r1",0)
            params.replace("three_prime_adapter", "TGGAATTCTCGGGTGCCAAGG")
        case 'nextflex':
            params.replace("clip_r1", 4)
            params.replace("three_prime_clip_r1", 4)
            params.replace("three_prime_adapter", "TGGAATTCTCGGGTGCCAAGG")
        case 'qiaseq':
            params.replace("clip_r1",0)
            params.replace("three_prime_clip_r1",0)
            params.replace("three_prime_adapter","AACTGTAGGCACCATCAAT")
        case 'cats':
            params.replace("clip_r1",3)
            params.replace("three_prime_clip_r1", 0)
            params.replace("three_prime_adapter", "AAAAAAAA")
        default: log.warn("No protocol specified, please ensure that you specified parameters for clipping/trimming or otherwise only auto-detection of adapters will be performed.")
    }
}
