// Import generic module functions
include { saveFiles; initOptions; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PARSE_FASTA_MIRNA {
    label 'process_medium'

    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::seqkit=2.0.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:2.0.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:2.0.0--h9ee0642_0"
    }

    input:
    path fasta

    output:
    path '*_igenome.fa' , emit: parsed_fasta

    script:
    def software = getSoftwareName(task.process)
    """
    # Uncompress FASTA reference files if necessary
    FASTA="$fasta"
    if [ \${FASTA: -3} == ".gz" ]; then
        gunzip -f \$FASTA
        FASTA=\${FASTA%%.gz}
    fi
    # Remove spaces from miRBase FASTA files
    # sed -i 's, ,_,g' \$FASTA
    sed '/^[^>]/s/[^AUGCaugc]/N/g' \$FASTA > fasta_parsed.fa
    sed -i 's/\s.*//' fasta_parsed.fa
    seqkit grep -r --pattern \".*${params.mirtrace_species}-.*\" fasta_parsed.fa > fasta_sps.fa
    seqkit seq --rna2dna fasta_sps.fa > fasta_igenome.fa

    #echo \$(seqkit --version 2>&1) | sed 's/^.*version //; s/Last.*\$//' > ${software}.version.txt
    """

}
