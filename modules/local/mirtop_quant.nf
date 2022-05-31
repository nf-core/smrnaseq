process MIRTOP_QUANT {
    label 'process_medium'

    conda (params.enable_conda ? 'mirtop=0.4.25 bioconda::samtools=1.15.1 conda-base::r-base=4.1.1 conda-base::r-data.table=1.14.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:63be019f50581cc5dfe4fc0f73ae50f2d4d661f7-0' :
        'quay.io/biocontainers/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:63be019f50581cc5dfe4fc0f73ae50f2d4d661f7-0' }"

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
    ${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

}
