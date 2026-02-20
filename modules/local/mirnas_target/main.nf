process MIRNA_TARGETS {
    label "process_medium"

    container 'docker.io/nfdata/multimir:v1.28.0.edger'

    input:
    path rdata
    tuple val (mirtrace_species), val (drug), val (disease) 

    output:
    path 'targets_*'                ,     emit: mirnas_targets
    path 'miRNA_targetDB_info.txt'  ,     emit: db_targets_info
    path 'miRNAs_targets.RData'     ,     emit: rdata         
    path('versions.yml')            ,     emit: versions
    path 'log.out'                  ,     emit: logfile
    
    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    
    """
    mirnas_targets.R $args ${rdata} ${mirtrace_species} ${drug} ${disease} > log.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        multiMiR: \$(Rscript -e "library(multiMiR); cat(as.character(packageVersion('multiMiR')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
    END_VERSIONS
    """

    stub:
    """
    touch miRNAs_targets.RData
    mkdir targets
    touch targets/predicted_DGE_summary.xlsx
    touch targets/predicted_DGE.xlsx
    touch targets/validated_DGE_summary.xlsx
    touch targets/validated_DGE.xlsx
    touch targets/validated_Downregulated_DGE.xlsx
    touch targets/validated_Upregulated_DGE.xlsx
    touch targets/miRNA_targetDB_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        multiMiR: \$(Rscript -e "library(multiMiR); cat(as.character(packageVersion('multiMiR')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
    END_VERSIONS
    """

}
