process RENAME_READS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b244720eef0bd28a41d4f26e33d3800e75d9fc87f080e81d42d5d676b4960dc/data' :
        'community.wave.seqera.io/library/gawk:5.3.0--180f75ae8b0ce739' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("renamed/*"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    // Without ext.prefix the reads keep their original file name
    def prefix = task.ext.prefix ?: ''
    """
    reads=( input/* )
    case "\${reads[0]}" in
        *.gz) extension=fastq.gz ;;
        *)    extension=fastq ;;
    esac
    name=\$(basename "\${reads[0]}")
    [ -n "${prefix}" ] && name="${prefix}.\${extension}"
    mkdir -p renamed
    ln -s "../\${reads[0]}" "renamed/\${name}"
    """

    stub:
    """
    mkdir -p renamed
    echo "" | gzip > renamed/${task.ext.prefix ?: meta.id}.fastq.gz
    """
}
