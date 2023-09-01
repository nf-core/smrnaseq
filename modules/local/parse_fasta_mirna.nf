process PARSE_FASTA_MIRNA {
    label 'process_medium'

    conda 'bioconda::seqkit=2.3.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0' :
        'biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(meta2), path(fasta)

    output:
    path '*_igenome.fa', emit: parsed_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filter_species = params.mirgenedb ? params.mirgenedb_species : params.mirtrace_species
    """
    # Uncompress FASTA reference files if necessary
    FASTA="$fasta"
    if [ \${FASTA: -3} == ".gz" ]; then
        gunzip -f \$FASTA
        FASTA=\${FASTA%%.gz}
    fi
    # Remove spaces from miRBase FASTA files
    # sed -i 's, ,_,g' \$FASTA
    sed '#^[^>]#s#[^AUGCaugc]#N#g' \$FASTA > \${FASTA}_parsed.fa
    # TODO perl -ane 's/[ybkmrsw]/N/ig;print;' \${FASTA}_parsed_tmp.fa > \${FASTA}_parsed.fa

    sed -i 's#\s.*##' \${FASTA}_parsed.fa
    seqkit grep -r --pattern \".*${filter_species}-.*\" \${FASTA}_parsed.fa > \${FASTA}_sps.fa
    seqkit seq --rna2dna \${FASTA}_sps.fa > \${FASTA}_igenome.fa

    cat <<-END_VERSIONS > versions.yml
    ${task.process}":
        seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
