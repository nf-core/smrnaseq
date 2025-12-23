process EDGER_QC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-edger_bioconductor-limma_r-base_r-data.table_pruned:8f1a97febf73b41c' :
        'community.wave.seqera.io/library/bioconductor-edger_bioconductor-limma_r-base_r-data.table_pruned:d592e095d0f6a5f0' }"

    input:
    path input_files

    output:
    path '*.{txt,pdf,csv}', emit: edger_files
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    edgeR_miRBase.r $input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
        edgeR: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        gplots: \$(Rscript -e "library(gplots); cat(as.character(packageVersion('gplots')))")
        methods: \$(Rscript -e "library(methods); cat(as.character(packageVersion('methods')))")
        statmod: \$(Rscript -e "library(statmod); cat(as.character(packageVersion('statmod')))")
    END_VERSIONS
    """

    stub:
    """
    touch "mature_unmapped_read_counts.txt"
    touch "hairpin_unmapped_read_counts.txt"
    touch "mature_counts.txt"
    touch "hairpin_counts.txt"
    touch "mature_logtpm.txt"
    touch "hairpin_logtpm.txt"
    touch "mature_logtpm.csv"
    touch "hairpin_logtpm.csv"
    touch "mature_normalized_CPM.txt.csv"
    touch "hairpin_normalized_CPM.txt.csv"
    touch "mature_CPM_heatmap.pdf"
    touch "hairpin_CPM_heatmap.pdf"
    touch "mature_edgeR_MDS_plot.pdf"
    touch "hairpin_edgeR_MDS_plot.pdf"
    touch "mature_edgeR_MDS_distance_matrix.txt"
    touch "hairpin_edgeR_MDS_distance_matrix.txt"
    touch "mature_edgeR_MDS_plot_coordinates.txt"
    touch "hairpin_edgeR_MDS_plot_coordinates.txt"
    touch "mature_log2CPM_sample_distances_heatmap.pdf"
    touch "hairpin_log2CPM_sample_distances_heatmap.pdf"
    touch "mature_log2CPM_sample_distances_dendrogram.pdf"
    touch "hairpin_log2CPM_sample_distances_dendrogram.pdf"
    touch "mature_log2CPM_sample_distances.txt"
    touch "hairpin_log2CPM_sample_distances.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
        edgeR: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')))")
        gplots: \$(Rscript -e "library(gplots); cat(as.character(packageVersion('gplots')))")
        methods: \$(Rscript -e "library(methods); cat(as.character(packageVersion('methods')))")
        statmod: \$(Rscript -e "library(statmod); cat(as.character(packageVersion('statmod')))")
    END_VERSIONS
    """
}
