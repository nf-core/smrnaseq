// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process EDGER {
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"edger", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::bioconductor-limma=3.50.0 bioconda::bioconductor-edger=3.36.0 conda-base::r-data.table=1.12.2 conda-base::r-gplot=3.1.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-limma:3.50.0--r41hd029910_0"
    } else {
        container "quay.io/biocontainers/bioconductor-limma:3.50.0--r41hd029910_0"
    }

    input:
    path input_files

    output:
    path '*.{txt,pdf,csv}' , emit: edger_files
    path "*version.txt" , emit: versions

    script:
    def software = getSoftwareName(task.process)
    """
    edgeR_miRBase.r $input_files
    cat <<-END_VERSIONS >> edgeR.version.txt
    \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    cat <<-END_VERSIONS >> limma.version.txt
    \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
    END_VERSIONS
    """

}
