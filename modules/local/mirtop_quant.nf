process MIRTOP_QUANT {
    label 'process_medium'

    conda 'mirtop=0.4.25 bioconda::samtools=1.15.1 conda-base::r-base=4.1.1 conda-base::r-data.table=1.14.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:63be019f50581cc5dfe4fc0f73ae50f2d4d661f7-0' :
        'biocontainers/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:63be019f50581cc5dfe4fc0f73ae50f2d4d661f7-0' }"

    input:
    path ("bams/*")
    path hairpin
    path gtf

    output:
    path "mirtop/mirtop.gff"        , emit: mirtop_gff
    path "mirtop/mirtop.tsv"        , emit: mirtop_table
    path "mirtop/mirtop_rawData.tsv", emit: mirtop_rawdata
    path "mirtop/stats/*"           , emit: logs
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filter_species = params.mirgenedb ? params.mirgenedb_species : params.mirtrace_species
    """
    #Cleanup the GTF if mirbase html form is broken
    GTF="$gtf"
    sed 's/&gt;/>/g' \$GTF | sed 's#<br>#\\n#g' | sed 's#</p>##g' | sed 's#<p>##g' | sed -e :a -e '/^\\n*\$/{\$d;N;};/\\n\$/ba' > \${GTF}_html_cleaned.gtf
    mirtop gff --hairpin $hairpin --gtf \${GTF}_html_cleaned.gtf -o mirtop --sps $filter_species ./bams/*
    mirtop counts --hairpin $hairpin --gtf \${GTF}_html_cleaned.gtf -o mirtop --sps $filter_species --add-extra --gff mirtop/mirtop.gff
    mirtop export --format isomir --hairpin $hairpin --gtf \${GTF}_html_cleaned.gtf --sps $filter_species -o mirtop mirtop/mirtop.gff
    mirtop stats mirtop/mirtop.gff --out mirtop/stats
    mv mirtop/stats/mirtop_stats.log mirtop/stats/full_mirtop_stats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

}
