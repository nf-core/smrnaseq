//
// This file holds several functions specific to the workflow/smrnaseq.nf in the nf-core/smrnaseq pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowSmrnaseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)
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
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
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
    public static String formatProtocol(params,log) {

        switch(params.protocol){
            case 'illumina':
                params.putIfAbsent("clip_r1", 1);
                params.putIfAbsent("three_prime_clip_r1",2);
                params.putIfAbsent("three_prime_adapter", "TGGAATTCTCGGGTGCCAAGG");
                break
            case 'nextflex':
                params.putIfAbsent("clip_r1", 4);
                params.putIfAbsent("three_prime_clip_r1", 4);
                params.putIfAbsent("three_prime_adapter", "TGGAATTCTCGGGTGCCAAGG");
                break
            case 'qiaseq':
                params.putIfAbsent("clip_r1",0);
                params.putIfAbsent("three_prime_clip_r1",0);
                params.putIfAbsent("three_prime_adapter","AACTGTAGGCACCATCAAT");
                break
            case 'cats':
                params.putIfAbsent("clip_r1",3);
                params.putIfAbsent("three_prime_clip_r1", 0);
                params.putIfAbsent("three_prime_adapter", "AAAAAAAA");
                break
            case 'custom':
                params.putIfAbsent("clip_r1", params.clip_r1)
                params.putIfAbsent("three_prime_clip_r1", params.three_prime_clip_r1)
            default:
                log.warn "Please make sure to specify all required clipping and trimming parameters, otherwise only adapter detection will be performed."
            }

            log.warn "Running with Protocol ${params.protocol}"
            log.warn "Therefore using Adapter: ${params.three_prime_adapter}"
            log.warn "Clipping ${params.clip_r1} bases from R1"
            log.warn "And clipping ${params.three_prime_clip_r1} bases from 3' end"
        }
}
