process ENRICHMENT {
    label "process_medium"

    //container 'docker.io/nfdata/bulk_rnaseq:v1.0.1'
    container 'docker.io/nfdata/clusterprofiler:v4.14.6'

    input:
    path rdata
    val mirtrace_species

    output: 
    path 'enrichment'           ,     emit: enrichr_results
    path 'orgdb_metadata.csv'   ,     emit: orgdb_metadata
    path 'enrichr.RData'        ,     emit: enrichr_RData
    path('versions.yml')        ,     emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    clusterprofiler.R "${rdata}" "${mirtrace_species}" > log.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        enrichR: \$(Rscript -e "library(enrichR); cat(as.character(packageVersion('enrichR')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
        enrichplot: \$(Rscript -e "library(enrichplot); cat(as.character(packageVersion('enrichplot')))")
        clusterProfiler: \$(Rscript -e "library(clusterProfiler); cat(as.character(packageVersion('clusterProfiler')))")
        AnnotationHub: \$(Rscript -e "library(AnnotationHub); cat(as.character(packageVersion('AnnotationHub')))")
        ggarchery: \$(Rscript -e "library(ggarchery); cat(as.character(packageVersion('ggarchery')))")

    END_VERSIONS
    """  
    
    stub:
    """
    mkdir enrichment
    mkdir all_targets
    touch orgdb_metadata.csv
    touch enrichment/all_targets/KEGG.xlsx
    touch enrichment/all_targets/GO_BiolProc.xlsx
    touch enrichment/all_targets/GO_CelCom.xlsx
    touch enrichment/all_targets/GO_MolFun.xlsx
    touch enrichment/all_targets/GO_BiolProc_barplot.png
    touch enrichment/all_targets/GO_CelCom_barplot.png
    touch enrichment/all_targets/GO_MolFun_barplot.png
    touch enrichment/all_targets/GO_BiolProc_dotplot.png
    touch enrichment/all_targets/GO_CelCom_dotplot.png
    touch enrichment/all_targets/GO_MolFun_dotplot.png
    touch enrichment/all_targets/GO_BiolProc_treeplot.png
    touch enrichment/all_targets/GO_CelCom_treeplot.png
    touch enrichment/all_targets/GO_MolFun_treeplot.png
    touch enrichment/all_targets/GO_BiolProc_emapplot.png
    touch enrichment/all_targets/GO_CelCom_emapplot.png
    touch enrichment/all_targets/GO_MolFun_emapplot.png
    touch enrichment/all_targets/KEGG_dotplot.png
    touch enrichment/all_targets/KEGG_dotplot.png

    mkdir downregulated
    touch enrichment/downregulated/KEGG.xlsx
    touch enrichment/downregulated/GO_BiolProc.xlsx
    touch enrichment/downregulated/GO_CelCom.xlsx
    touch enrichment/downregulated/GO_MolFun.xlsx
    touch enrichment/downregulated/GO_BiolProc_barplot.png
    touch enrichment/downregulated/GO_CelCom_barplot.png
    touch enrichment/downregulated/GO_MolFun_barplot.png
    touch enrichment/downregulated/GO_BiolProc_dotplot.png
    touch enrichment/downregulated/GO_CelCom_dotplot.png
    touch enrichment/downregulated/GO_MolFun_dotplot.png
    touch enrichment/downregulated/GO_BiolProc_treeplot.png
    touch enrichment/downregulated/GO_CelCom_treeplot.png
    touch enrichment/downregulated/GO_MolFun_treeplot.png
    touch enrichment/downregulated/GO_BiolProc_emapplot.png
    touch enrichment/downregulated/GO_CelCom_emapplot.png
    touch enrichment/downregulated/GO_MolFun_emapplot.png
    touch enrichment/downregulated/KEGG_dotplot.png
    touch enrichment/downregulated/KEGG_dotplot.png

    mkdir upregulated
    touch enrichment/upregulated/KEGG.xlsx
    touch enrichment/upregulated/GO_BiolProc.xlsx
    touch enrichment/upregulated/GO_CelCom.xlsx
    touch enrichment/upregulated/GO_MolFun.xlsx
    touch enrichment/upregulated/GO_BiolProc_barplot.png
    touch enrichment/upregulated/GO_CelCom_barplot.png
    touch enrichment/upregulated/GO_MolFun_barplot.png
    touch enrichment/upregulated/GO_BiolProc_dotplot.png
    touch enrichment/upregulated/GO_CelCom_dotplot.png
    touch enrichment/upregulated/GO_MolFun_dotplot.png
    touch enrichment/upregulated/GO_BiolProc_treeplot.png
    touch enrichment/upregulated/GO_CelCom_treeplot.png
    touch enrichment/upregulated/GO_MolFun_treeplot.png
    touch enrichment/upregulated/GO_BiolProc_emapplot.png
    touch enrichment/upregulated/GO_CelCom_emapplot.png
    touch enrichment/upregulated/GO_MolFun_emapplot.png
    touch enrichment/upregulated/KEGG_dotplot.png
    touch enrichment/upregulated/KEGG_dotplot.png


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        enrichR: \$(Rscript -e "library(enrichR); cat(as.character(packageVersion('enrichR')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        openxlsx: \$(Rscript -e "library(openxlsx); cat(as.character(packageVersion('openxlsx')))")
        enrichplot: \$(Rscript -e "library(enrichplot); cat(as.character(packageVersion('enrichplot')))")
        clusterProfiler: \$(Rscript -e "library(clusterProfiler); cat(as.character(packageVersion('clusterProfiler')))")
        AnnotationHub: \$(Rscript -e "library(AnnotationHub); cat(as.character(packageVersion('AnnotationHub')))")
        ggarchery: \$(Rscript -e "library(ggarchery); cat(as.character(packageVersion('ggarchery')))")
    END_VERSIONS
    """
  
}