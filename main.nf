#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
          N G I    S M A L L    R N A - S E Q    B E S T    P R A C T I C E
========================================================================================
 Small-RNA-Seq Best Practice Analysis Pipeline. Started May 2016.
 #### Homepage / Documentation
 https://github.com/SciLifeLab/NGI-smRNAseq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
 Chuan Wang <chuan.wang@scilifelab.se>
 Rickard Hammarén <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   FastQC for raw sequencing reads quility control
 - 2:   Trim Galore! for adapter trimming
 - 3.1: Bowtie 1 alignment against miRBase mature miRNA
 - 3.2: Post-alignment processing of miRBase mature miRNA counts
 - 3.3: edgeR analysis on miRBase mature miRNA counts
        - TMM normalization and a table of top expression mature miRNA
        - MDS plot clustering samples
        - Heatmap of sample similarities
 - 4.1: Bowtie 1 alignment against miRBase hairpin for the unaligned reads in step 3
 - 4.2: Post-alignment processing of miRBase hairpin counts
 - 4.3: edgeR analysis on miRBase hairpin counts
        - TMM normalization and a table of top expression hairpin
        - MDS plot clustering samples
        - Heatmap of sample similarities
 - 5.1: Bowtie 2 alignment against host reference genome
 - 5.2: Post-alignment processing of Bowtie 2
 - 6:   NGI-Visualization of Bowtie 2 alignment statistics
 - 7:   MultiQC
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     NGI-smRNAseq : smRNA-Seq Best Practice v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run SciLifeLab/NGI-smRNAseq --reads '*.fastq.gz' --genome GRCh37

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
                                    NOTE! Paired-end data is NOT supported by this pipeline! For paired-end data, use Read 1 only.
      --genome                      Name of iGenomes reference
                                    NOTE! With the option --genome 'ALL', the entire dataset of mature miRNAs and hairpins
                                    in miRBase will be used as reference regardless of species. Meanwhile the alignment against
                                    host reference genome will be skipped.

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --saveReference               Save the generated reference files the the Results directory.

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --rlocation                   Location to save R-libraries used in the pipeline. Default value is ~/R/nxtflow_libs/
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.5

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

// Configurable variables
params.name = false
params.project = false
params.genome = false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bt2index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.mature = params.genome ? params.genomes[ params.genome ].mature ?: false : false
params.hairpin = params.genome ? params.genomes[ params.genome ].hairpin ?: false : false
params.saveReference = false
params.reads = "data/*.fastq.gz"
params.outdir = './results'
params.email = false
params.plaintext_email = false

// R library locations
params.rlocation = false
if (params.rlocation){
    nxtflow_libs = file(params.rlocation)
    nxtflow_libs.mkdirs()
}

// Validate inputs
if( !params.mature || !params.hairpin ){
    exit 1, "Missing mature / hairpin reference indexes! Is --genome specified?"
}
if( params.mature ){
    mature = file(params.mature)
    if( !mature.exists() ) exit 1, "Mature file not found: ${params.mature}"
}
if( params.hairpin ){
    hairpin = file(params.hairpin)
    if( !hairpin.exists() ) exit 1, "Hairpin file not found: ${params.hairpin}"
}
if( params.gtf ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}
if( params.bt2index ){
    bt2_index = file("${params.bt2index}.1.bt2")
    bt2_indices = Channel.fromPath( "${params.bt2index}*" ).toList()
    if( !bt2_index.exists() ) exit 1, "Reference genome Bowtie 2 not found: ${params.bt2index}"
}
if( !params.gtf || !params.bt2index) {
    log.info "No GTF / Bowtie2 index supplied - host reference genome analysis will be skipped."
}
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { raw_reads_fastqc; raw_reads_trimgalore }

// Header log info
log.info "==========================================="
log.info " NGI-smRNAseq : Small RNA-Seq Best Practice v${version}"
log.info "==========================================="
def summary = [:]
summary['Run Name']            = custom_runName ?: workflow.runName
summary['Reads']               = params.reads
summary['Genome']              = params.genome
summary['miRBase mature']      = params.mature
summary['miRBase hairpin']     = params.hairpin
if(params.bt2index)            summary['Bowtie2 Index'] = params.bt2index
if(params.gtf)                 summary['GTF Annotation'] = params.gtf
summary['Save Reference']      = params.saveReference ? 'Yes' : 'No'
summary['Output dir']          = params.outdir
summary['Working dir']         = workflow.workDir
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['R libraries']         = params.rlocation
summary['Script dir']          = workflow.projectDir
summary['Config Profile'] = (workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile)
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================="


/*
 * PREPROCESSING - Build Bowtie index for mature and hairpin
 */
process makeBowtieIndex {

    publishDir path: { params.saveReference ? "${params.outdir}/bowtie/reference" : params.outdir },
               saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file mature from mature
    file hairpin from hairpin

    output:
    file 'mature_idx.*' into mature_index
    file 'hairpin_idx.*' into hairpin_index

    script:
    """
    fasta_formatter -w 0 -i $mature -o mature_igenome.fa
    fasta_nucleotide_changer -d -i mature_igenome.fa -o mature_idx.fa
    bowtie-build mature_idx.fa mature_idx
    fasta_formatter -w 0 -i $hairpin -o hairpin_igenome.fa
    fasta_nucleotide_changer -d -i hairpin_igenome.fa -o hairpin_idx.fa
    bowtie-build hairpin_idx.fa hairpin_idx
    """
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$reads"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file reads from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -q $reads
    fastqc --version
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$reads"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    file reads from raw_reads_trimgalore

    output:
    file '*.gz' into trimmed_reads_bowtie, trimmed_reads_bowtie2, trimmed_reads_insertsize
    file '*trimming_report.txt' into trimgalore_results, trimgalore_logs
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
    """
    trim_galore --small_rna --gzip $reads --fastqc
    """
}


/*
 * STEP 2.1 - Insertsize
 */

process insertsize {
    tag "$name"
    publishDir "${params.outdir}/trim_galore/insertsize", mode: 'copy'

    input:
    set val(name), file(reads) from trimmed_reads_insertsize

    output:
    file '*.insertsize' into insertsize_results

    script:
    """
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' <(zcat \$input) >\${prefix}.insertsize
    """
}


/*
 * STEP 3 - Bowtie miRBase mature miRNA
 */
process bowtie_miRBase_mature {
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_mature", mode: 'copy', pattern: '*.mature_unmapped.fq.gz'

    input:
    file reads from trimmed_reads_bowtie
    file index from mature_index

    output:
    file '*.mature.bam' into miRBase_mature_bam
    file '*.mature_unmapped.fq.gz' into mature_unmapped_reads
    file '.command.log' into bowtie_log, bowtie_mature_alignment

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bowtie \\
        $index_base \\
        -q <(zcat $reads) \\
        -p 2 \\
        -t \\
        -n 0 \\
        -l 15 \\
        -e 99999 \\
        -k 10 \\
        --best \\
        --chunkmbs 2048 \\
        --un ${prefix}.mature_unmapped.fq \\
        -S \\
        | samtools view -bS - > ${prefix}.mature.bam
    bowtie --version

    gzip ${prefix}.mature_unmapped.fq
    """
}

/*
 * STEP 4 - Bowtie against miRBase hairpin
 */
process bowtie_miRBase_hairpin {
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_hairpin", mode: 'copy', pattern: '*.hairpin_unmapped.fq.gz'

    input:
    file reads from mature_unmapped_reads
    file index from hairpin_index

    output:
    file '*.hairpin.bam' into miRBase_hairpin_bam
    file '*.hairpin_unmapped.fq.gz' into hairpin_unmapped_reads
    file '.command.log' into bowtie_hairpin_alignment

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - '.mature_unmapped.fq.gz'
    """
    bowtie \\
        $index_base \\
        -p 2 \\
        -t \\
        -n 1 \\
        -l 15 \\
        -e 99999 \\
        -k 10 \\
        --best \\
        --chunkmbs 2048 \\
        -q <(zcat $reads) \\
        --un ${prefix}.hairpin_unmapped.fq \\
        -S \\
        | samtools view -bS - > ${prefix}.hairpin.bam

    gzip ${prefix}.hairpin_unmapped.fq
    """
}


/*
 * STEP 5 - Post-alignment processing for miRBase mature and hairpin
 */
def wrap_mature_and_hairpin = { file ->
    if ( file.contains("mature") ) return "miRBase_mature/$file"
    if ( file.contains("hairpin") ) return "miRBase_hairpin/$file"
}

process miRBasePostAlignment {
    tag "$input"
    publishDir "${params.outdir}/bowtie", mode: 'copy', saveAs: wrap_mature_and_hairpin

    input:
    file input from miRBase_mature_bam.mix(miRBase_hairpin_bam)

    output:
    file "${input.baseName}.count" into miRBase_counts
    file "${input.baseName}.sorted.bam" into miRBase_bam
    file "${input.baseName}.sorted.bam.bai" into miRBase_bai

    script:
    """
    samtools sort ${input.baseName}.bam -o ${input.baseName}.sorted.bam
    samtools index ${input.baseName}.sorted.bam
    samtools idxstats ${input.baseName}.sorted.bam > ${input.baseName}.count
    """
}


/*
 * STEP 6 - edgeR miRBase feature counts processing
 */
process edgeR_miRBase {
    publishDir "${params.outdir}/edgeR", mode: 'copy', saveAs: wrap_mature_and_hairpin

    input:
    file input_files from miRBase_counts.toSortedList()

    output:
    file '*.{txt,pdf}' into edgeR_miRBase_results

    script:
    """
    edgeR_miRBase.r $params.rlocation $input_files
    """
}


/*
 * STEP 7.1 and 7.2 IF A GENOME SPECIFIED ONLY!
 */
if( params.gtf && params.bt2index) {

    /*
     * STEP 7.1 - Bowtie 2 against reference genome
     */
    process bowtie2 {
        tag "$reads"
        publishDir "${params.outdir}/bowtie2", mode: 'copy'

        input:
        file reads from trimmed_reads_bowtie2
        file index from bt2_index
        file bt2_indices

        output:
        file '*.bowtie2.bam' into bowtie2_bam, bowtie2_bam_for_unmapped
        file '.command.log' into bowtie2_log, bowtie2_alignment

        script:
        index_base = index.toString() - '.1.bt2'
        prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        bowtie2 \\
            -x $index_base \\
            -U $reads \\
            -k 10 \\
            --very-sensitive \\
            -p 8 \\
            -t \\
            | samtools view -bT $index_base - > ${prefix}.bowtie2.bam
        bowtie2 --version
        """
    }

    /*
     * STEP 7.2 - Bowtie 2 Statistics about unmapped reads against ref genome
     */

    process bowtie2_unmapped {

        publishDir "${params.outdir}/bowtie2/unmapped", mode: 'copy'

        input:
        file input_files from bowtie2_bam_for_unmapped.toSortedList()

        output:
        file 'unmapped_refgenome.txt' into bowtie2_unmapped

        script:
        """
        for i in $input_files
        do
          printf "\${i}\t"
          samtools view -c -f0x4 \${i}
        done > unmapped_refgenome.txt
        """
    }


    /*
     * STEP 7.3 - NGI-Visualizations of Bowtie 2 alignment statistics
     */
    process ngi_visualizations {
        tag "$bowtie2_bam"
        publishDir "${params.outdir}/bowtie2/ngi_visualizations", mode: 'copy'

        input:
        file gtf from gtf
        file bowtie2_bam

        output:
        file '*.{png,pdf}' into bowtie2_ngi_visualizations

        script:
        // Note! ngi_visualizations needs to be installed!
        // See https://github.com/NationalGenomicsInfrastructure/ngi_visualizations
        """
        #!/usr/bin/env python
        from ngi_visualizations.biotypes import count_biotypes
        count_biotypes.main('$gtf','$bowtie2_bam')
        """
    }

}

/*
 * Parse software version numbers
 */
software_versions = [
  'FastQC': null, 'Trim Galore!': null, 'Bowtie': null, 'Nextflow': "v$workflow.nextflow.version"
]
if( params.gtf && params.bt2index ) software_versions['Bowtie 2'] = null

process get_software_versions {
    cache false
    executor 'local'

    input:
    val fastqc from fastqc_stdout.collect()
    val trim_galore from trimgalore_logs.collect()
    val bowtie from bowtie_log.collect()
    val bowtie2 from bowtie2_log.collect()

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    exec:
    software_versions['FastQC'] = fastqc[0].getText().find(/FastQC v(\S+)/) { match, version -> "v$version" }
    software_versions['Trim Galore!'] = trim_galore[0].getText().find(/Trim Galore version: (\S+)/) {match, version -> "v$version"}
    software_versions['Bowtie'] = bowtie[0].getText().find(/bowtie-align version (\S+)/) { match, version -> "v$version" }
    if( software_versions.containsKey('Bowtie') ) software_versions['Bowtie 2'] = bowtie2[0].getText().find(/bowtie2-align-s version (\S+)/) { match, version -> "v$version" }

    def sw_yaml_file = task.workDir.resolve('software_versions_mqc.yaml')
    sw_yaml_file.text  = """
    id: 'ngi-smrnaseq'
    section_name: 'NGI-smRNAseq Software Versions'
    section_href: 'https://github.com/SciLifeLab/NGI-smRNAseq'
    plot_type: 'html'
    description: 'are collected at run time from the software output.'
    data: |
        <dl class=\"dl-horizontal\">
${software_versions.collect{ k,v -> "            <dt>$k</dt><dd>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</dd>" }.join("\n")}
        </dl>
    """.stripIndent()
}

/*
 * STEP 8 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trim_galore/*') from trimgalore_results.flatten().toList()
    file ('trim_galore/*') from trimgalore_fastqc_reports.flatten().toList()
    file ('bowtie/miRBase_mature/*') from bowtie_mature_alignment.flatten().toList()
    file ('bowtie/miRBase_hairpin/*') from bowtie_haripin_alignment.flatten().toList()
    file ('edgeR/*') from edgeR_miRBase_results.flatten().toList()
    if( params.gtf && params.bt2index ) file ('bowtie2/*') from bowtie2_alignment.flatten().toList()
    file ('software_versions/*') from software_versions_yaml

    output:
    file '*multiqc_report.html' into multiqc_html
    file '*multiqc_data' into multiqc_data

    script:
    """
    multiqc -f .
    """
}
multiqc_stderr.subscribe { stderr ->
  software_versions['MultiQC'] = stderr.getText().find(/This is MultiQC v(\S+)/) { match, version -> "v$version" }
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[NGI-smRNAseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[NGI-smRNAseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['software_versions'] = software_versions
    email_fields['software_versions']['Nextflow Build'] = workflow.nextflow.build
    email_fields['software_versions']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[NGI-smRNAseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[NGI-smRNAseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    ngismrnaseqlogo = new File("$baseDir/assets/NGI-smRNAseq_logo.png").bytes.encodeBase64().toString()
    scilifelablogo = new File("$baseDir/assets/SciLifeLab_logo.png").bytes.encodeBase64().toString()
    ngilogo = new File("$baseDir/assets/NGI_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:ngismrnaseqlogo/, "data:image/png;base64,$ngismrnaseqlogo")
    email_html = email_html.replaceAll(~/cid:scilifelablogo/, "data:image/png;base64,$scilifelablogo")
    email_html = email_html.replaceAll(~/cid:ngilogo/, "data:image/png;base64,$ngilogo")

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NGI-smRNAseq] Pipeline Complete"

}
