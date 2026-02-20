process MIRNAS_INTEREST {
    container 'docker.io/nfdata/multimir:v1.28.0.edger'
    tag "process_single"

    input:
    path rdata
    path mirnaslist
    val mirtrace_species

    // Define the output files
    output:
    path '*.{pdf,png}'            ,     emit: mirnas_interest_plots
    path '*.xlsx'                 ,     emit: mirnas_interest
    path('versions.yml')          ,     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mirnas_interest.R "${rdata}" "${mirnaslist}" "${mirtrace_species}" > log.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        multiMiR: \$(Rscript -e "library(multiMiR); cat(as.character(packageVersion('multiMiR')))")
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
        RColorBrewer: \$(Rscript -e "library(RColorBrewer); cat(as.character(packageVersion('RColorBrewer')))")
        gplots: \$(Rscript -e "library(gplots); cat(as.character(packageVersion('gplots')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir mirnas_interest
    touch mirnas_interest/genes_interest_summary.xlsx
    touch mirnas_interest/genes_interest.xlsx
    touch mirnas_interest/mirnas_interest_summary.xlsx
    touch mirnas_interest/mirnas_interest.xlsx
    touch mirnas_interest/heatmap.pdf
    touch mirnas_interest/heatmap.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        multiMiR: \$(Rscript -e "library(multiMiR); cat(as.character(packageVersion('multiMiR')))")
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
        RColorBrewer: \$(Rscript -e "library(RColorBrewer); cat(as.character(packageVersion('RColorBrewer')))")
        gplots: \$(Rscript -e "library(gplots); cat(as.character(packageVersion('gplots')))")
    END_VERSIONS
    """
}
