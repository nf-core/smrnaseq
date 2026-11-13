nextflow.enable.types = true

process EDGER_QC {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7b/7b75cd43418d68682b9872aafa0e3ef683cf609af2a70cf33f27a323e2d718c9/data' :
        'community.wave.seqera.io/library/bioconductor-edger_bioconductor-limma_r-base_r-data.table_pruned:d592e095d0f6a5f0' }"

    input:
    input_files: Set<Path>

    output:
    edger_files = files('*.{txt,pdf,csv}')

    topic:
    tuple(task.process, 'r-base', eval("Rscript -e 'cat(as.character(getRversion()))'")) >> 'versions'
    tuple(task.process, 'limma', eval("Rscript -e \"cat(as.character(packageVersion('limma')))\"")) >> 'versions'
    tuple(task.process, 'edgeR', eval("Rscript -e \"cat(as.character(packageVersion('edgeR')))\"")) >> 'versions'
    tuple(task.process, 'data.table', eval("Rscript -e \"cat(as.character(packageVersion('data.table')))\"")) >> 'versions'
    tuple(task.process, 'gplots', eval("Rscript -e \"cat(as.character(packageVersion('gplots')))\"")) >> 'versions'
    tuple(task.process, 'methods', eval("Rscript -e \"cat(as.character(packageVersion('methods')))\"")) >> 'versions'
    tuple(task.process, 'statmod', eval("Rscript -e \"cat(as.character(packageVersion('statmod')))\"")) >> 'versions'

    script:
    """
    edgeR_miRBase.r ${input_files.join(" ")}
    """

    stub:
    """
    touch "mature_unmapped_read_counts.txt"
    touch "hairpin_unmapped_read_counts.txt"
    touch "mature_counts.csv"
    touch "hairpin_counts.csv"
    touch "mature_logtpm.txt"
    touch "hairpin_logtpm.txt"
    touch "mature_logtpm.csv"
    touch "hairpin_logtpm.csv"
    touch "mature_normalized_CPM.txt"
    touch "hairpin_normalized_CPM.txt"
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
    """
}
