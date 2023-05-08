process EDGER_QC {
    label 'process_medium'

    conda 'bioconda::bioconductor-limma=3.50.0 bioconda::bioconductor-edger=3.36.0 conda-forge::r-data.table=1.14.2 conda-forge::r-gplots=3.1.1 conda-forge::r-statmod=1.4.36'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0' :
        'biocontainers/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0' }"

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

}
