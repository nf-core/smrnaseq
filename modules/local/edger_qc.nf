// Import generic module functions
include { saveFiles; initOptions; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process EDGER_QC {
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"edger", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::bioconductor-limma=3.50.0 bioconda::bioconductor-edger=3.36.0 conda-forge::r-data.table=1.14.2 conda-forge::r-gplots=3.1.1 conda-forge::r-statmod=1.4.36' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0"
    }

    input:
    path input_files

    output:
    path '*.{txt,pdf,csv}', emit: edger_files
    path "versions.yml"   , emit: versions

    script:
    """
    edgeR_miRBase.r $input_files

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
    END_VERSIONS
    """

}
