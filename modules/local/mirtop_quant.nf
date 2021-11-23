// Import generic module functions
include { saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

process MIRTOP_QUANT {
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:".", meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::mirtop=0.4.24 bioconda::samtools=1.13 conda-base::r-base=4.0.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:c8b426b3677a47f58bbbc443c0cf67bf8d53a97b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:c8b426b3677a47f58bbbc443c0cf67bf8d53a97b-0"
    }

    input:
    path ("bams/*")
    path hairpin
    path gtf

    output:
    path "mirtop/mirtop.gff"
    path "mirtop/mirtop.tsv"        , emit: mirtop_table
    path "mirtop/mirtop_rawData.tsv"
    path "mirtop/stats/*"           , emit: logs
    path "versions.yml"             , emit: versions

    script:
    """
    mirtop gff --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species ./bams/*
    mirtop counts --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species --add-extra --gff mirtop/mirtop.gff
    mirtop export --format isomir --hairpin $hairpin --gtf $gtf --sps $params.mirtrace_species -o mirtop mirtop/mirtop.gff
    mirtop stats mirtop/mirtop.gff --out mirtop/stats
    mv mirtop/stats/mirtop_stats.log mirtop/stats/full_mirtop_stats.log

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

}
