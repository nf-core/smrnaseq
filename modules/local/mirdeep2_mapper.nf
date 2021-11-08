// Import generic module functions
include { saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

def VERSION = '2.0.1'

process MIRDEEP2_MAPPER {
    label 'process_medium'
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::mirdeep2=2.0.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.3--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/mirdeep2:2.0.1.3--hdfd78af_1"
    }

    when:
    !params.skip_mirdeep  // TODO ? I think it would be better to have this logic outside the module

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple path('*_collapsed.fa'), path('*reads_vs_refdb.arf'), emit: mirdeep2_inputs
    path "versions.yml"                                      , emit: versions


    script:
    def index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    """
    mapper.pl \\
    $reads \\
    -e \\
    -h \\
    -i \\
    -j \\
    -m \\
    -p $index_base \\
    -s ${meta.id}_collapsed.fa \\
    -t ${meta.id}_reads_vs_refdb.arf \\
    -o 4

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo "$VERSION")
    END_VERSIONS
    """
}
