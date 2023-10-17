process BLAT_MIRNA {
    tag "$fasta"
    label 'process_medium'

    conda 'bioconda::blat=36'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blat:36--0' :
        'biocontainers/blat:36--0' }"

    input:
    val db_type
    path mirna
    path contaminants


    output:
    path 'filtered.fa'  , emit: filtered_set
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if ( db_type == "cdna" )
        """
        echo $db_type
        awk '/^>/ { x=index(\$6, "transcript_biotype:miRNA") } { if(!x) print }' $contaminants > subset.fa
        blat -out=blast8 $mirna subset.fa /dev/stdout | awk 'BEGIN{FS="\t"}{if(\$11 < 1e-5)print \$1;}' | uniq > mirnahit.txt
        awk 'BEGIN { while((getline<"mirnahit.txt")>0) l[">"\$1]=1 } /^>/ {x = l[\$1]} {if(!x) print }' subset.fa  > filtered.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blat: \$(echo \$(blat) | grep Standalone | awk '{ if (match(\$0,/[0-9]*[0-9]/,m)) print m[0] }')
        END_VERSIONS
        """

    else if ( db_type == "ncrna" )
        """
        echo $db_type
        awk '/^>/ { x=(index(\$6, "transcript_biotype:rRNA") || index(\$6, "transcript_biotype:miRNA")) } { if(!x) print }' $contaminants > subset.fa
        blat -out=blast8 $mirna subset.fa /dev/stdout | awk 'BEGIN{FS="\t"}{if(\$11 < 1e-5)print \$1;}' | uniq > mirnahit.txt
        awk 'BEGIN { while((getline<"mirnahit.txt")>0) l[">"\$1]=1 } /^>/ {x = l[\$1]} {if(!x) print }' subset.fa  > filtered.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blat: \$(echo \$(blat) | grep Standalone | awk '{ if (match(\$0,/[0-9]*[0-9]/,m)) print m[0] }')
        END_VERSIONS
        """

    else
        """
        echo $db_type
        blat -out=blast8 $mirna $contaminants /dev/stdout | awk 'BEGIN{FS="\t"}{if(\$11 < 1e-5)print \$1;}' | uniq > mirnahit.txt
        awk 'BEGIN { while((getline<"mirnahit.txt")>0) l[">"\$1]=1 } /^>/ {x = l[\$1]} {if(!x) print }' $contaminants  > filtered.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blat: \$(echo \$(blat) | grep Standalone | awk '{ if (match(\$0,/[0-9]*[0-9]/,m)) print m[0] }')
        END_VERSIONS
        """

}
