nextflow.enable.types = true

process PARSE_FASTA_MIRNA {
    label 'process_medium'

    conda 'bioconda::seqkit=2.6.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/913600c87cb9a19229a5718887849281df9b53e49121913574587a3368dc8dfc/data' :
        'community.wave.seqera.io/library/seqkit:2.6.1--49efc1ecf715e29f' }"

    input:
    tuple(meta2: Map, fasta: Path)
    filter_species: String

    output:
    parsed_fasta = tuple(meta2, file('*_igenome.fa'))

    topic:
    tuple(task.process, 'seqkit', eval("seqkit version | sed 's/^seqkit v//'")) >> 'versions'

    script:
    """
    # Uncompress FASTA reference files if necessary
    FASTA="$fasta"
    if [ \${FASTA: -3} == ".gz" ]; then
        gunzip -f \$FASTA
        FASTA=\${FASTA%%.gz}
    fi
    sed 's/&gt;/>/g' \$FASTA | sed 's#<br>#\\n#g' | sed 's#</p>##g' | sed 's#<p>##g' | sed -e :a -e '/^\\n*\$/{\$d;N;};/\\n\$/ba' > \${FASTA}_html_cleaned.fa
    # Remove spaces from miRBase FASTA files
    sed '#^[^>]#s#[^AUGCaugc]#N#g' \${FASTA}_html_cleaned.fa > \${FASTA}_parsed.fa

    sed -i 's#\s.*##' \${FASTA}_parsed.fa
    seqkit grep -r --pattern \".*${filter_species}-.*\" \${FASTA}_parsed.fa > \${FASTA}_sps.fa
    seqkit seq --rna2dna \${FASTA}_sps.fa > \${FASTA}_igenome.fa
    """

    stub:
    // Same output name as the script: a .gz suffix is dropped before adding _igenome.fa
    def fasta_name = fasta.name.endsWith('.gz') ? fasta.baseName : fasta.name
    """
    touch ${fasta_name}_igenome.fa
    """
}
