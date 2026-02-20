process DGE_EDGER {
    label "process_single"

    container 'docker.io/nfdata/edger:v4.4.2'

    input:
    path counts
    path metadata
    val contrasts
    val lfc
    val fdr

    output:
    path "DE_*"             , emit: edger_files
    path("DGE_miRNAs.RData"), emit: rdata
    path("*.{pdf,png}")      , emit: plots
    path('versions.yml')    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    """
    dge_edger.R $args "${counts}" "${metadata}" "${contrasts}" ${lfc} ${fdr} > log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
        edgeR: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
        RColorBrewer: \$(Rscript -e "library(RColorBrewer); cat(as.character(packageVersion('RColorBrewer')))")
        gplots: \$(Rscript -e "library(gplots); cat(as.character(packageVersion('gplots')))")
        EnhancedVolcano: \$(Rscript -e "library(EnhancedVolcano); cat(as.character(packageVersion('EnhancedVolcano')))")
    END_VERSIONS
    """

    stub:
    """
    touch DGE_miRNAs.RData
    mkdir dge_contrasts
    touch dge_contrasts/DGE_analysis_summary.txt
    touch dge_contrasts/DGE.xlsx
    touch dge_contrasts/normalized_CPM.xlsx
    mkdir dge_contrasts/plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
        edgeR: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
        RColorBrewer: \$(Rscript -e "library(RColorBrewer); cat(as.character(packageVersion('RColorBrewer')))")
        gplots: \$(Rscript -e "library(gplots); cat(as.character(packageVersion('gplots')))")
        EnhancedVolcano: \$(Rscript -e "library(EnhancedVolcano); cat(as.character(packageVersion('EnhancedVolcano')))")
    END_VERSIONS
    """
}
