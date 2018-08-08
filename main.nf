#!/usr/bin/env nextflow

/*
========================================================================================
               S M A L L    R N A - S E Q    B E S T    P R A C T I C E
========================================================================================
 Small-RNA-Seq Best Practice Analysis Pipeline. Started May 2016.
 #### Homepage / Documentation
 https://github.com/nf-core/smrnaseq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
 Chuan Wang <chuan.wang@scilifelab.se>
 Rickard Hammarén <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   FastQC for raw sequencing reads quality control
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
 - 5.1: Bowtie 1 alignment against host reference genome
 - 5.2: Post-alignment processing of Bowtie 1 alignment against host reference genome
 - 6:   NGI-Visualization of Bowtie alignment statistics
 - 7:   MultiQC
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     nf-core/smrnaseq : smRNA-Seq Best Practice v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/smrnaseq --reads '*.fastq.gz' --genome GRCh37

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
                                    NOTE! Paired-end data is NOT supported by this pipeline! For paired-end data, use Read 1 only
      --genome                      Name of iGenomes reference

    References
      --saveReference               Save the generated reference files the the Results directory
      --mature                      Path to the FASTA file of mature miRNAs
      --hairpin                     Path to the FASTA file of miRNA precursors
      --bt_index                    Path to the bowtie 1 index files of the host reference genome

    Trimming options
      --length [int]                Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18
      --clip_R1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1
      --three_prime_clip_R1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
      --seqCenter                   Text about sequencing center which will be added in the header of output bam files (Note that no blank is allowed!)
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Pipeline options
params.name = false
params.project = false
params.genome = false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bt_index = params.genome ? params.genomes[ params.genome ].bowtie ?: false : false
params.bt_indices = null
params.mature = params.genome ? params.genomes[ params.genome ].mature ?: false : false
params.hairpin = params.genome ? params.genomes[ params.genome ].hairpin ?: false : false
params.saveReference = false
params.reads = "data/*.fastq.gz"
params.readPaths = null
params.email = false
params.plaintext_email = false
params.seqCenter = false

// Custom trimming options
params.length = 18
params.clip_R1 = 0
params.three_prime_clip_R1 = 0

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
if( params.bt_index ){
    bt_index = file("${params.bt_index}.fa")
    bt_indices = Channel.fromPath( "${params.bt_index}*.ebwt" ).toList()
    if( !bt_index.exists() ) exit 1, "Reference genome for Bowtie 1 not found: ${params.bt_index}"
} else if( params.bt_indices ){
    bt_indices = Channel.from(params.readPaths).map{ file(it) }.toList()
}
if( !params.gtf || !params.bt_index) {
    log.info "No GTF / Bowtie 1 index supplied - host reference genome analysis will be skipped."
}
multiqc_config = file(params.multiqc_config)

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */
if(params.readPaths){
    Channel
        .from(params.readPaths)
        .map { file(it) }
        .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        .into { raw_reads_fastqc; raw_reads_trimgalore }
} else {
    Channel
        .fromPath( params.reads )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
        .into { raw_reads_fastqc; raw_reads_trimgalore }
}


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/smrnaseq : Small RNA-Seq Best Practice v${params.version}
======================================================="""
def summary = [:]
summary['Run Name']            = custom_runName ?: workflow.runName
summary['Reads']               = params.reads
summary['Genome']              = params.genome
summary['Trim min length']     = params.length
summary["Trim 5' R1"]          = params.clip_R1
summary["Trim 3' R1"]          = params.three_prime_clip_R1
summary['miRBase mature']      = params.mature
summary['miRBase hairpin']     = params.hairpin
if(params.bt_index)            summary['Bowtie Index for Ref'] = params.bt_index
if(params.gtf)                 summary['GTF Annotation'] = params.gtf
summary['Save Reference']      = params.saveReference ? 'Yes' : 'No'
summary['Output dir']          = params.outdir
summary['Working dir']         = workflow.workDir
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Script dir']          = workflow.projectDir
summary['Config Profile'] = (workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile)
if(params.project) summary['UPPMAX Project'] = params.project
if(params.seqCenter) summary['Seq Center'] = params.seqCenter
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}
// Show a big error message if we're running on the base config and an uppmax cluster
if( workflow.profile == 'standard'){
    if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
        log.error "====================================================\n" +
                  "  WARNING! You are running with the default 'standard'\n" +
                  "  pipeline config profile, which runs on the head node\n" +
                  "  and assumes all software is on the PATH.\n" +
                  "  ALL JOBS ARE RUNNING LOCALLY and stuff will probably break.\n" +
                  "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                  "============================================================"
    }
}


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

    script:
    """
    fastqc -q $reads
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
    file '*.gz' into trimmed_reads_bowtie, trimmed_reads_bowtie_ref, trimmed_reads_insertsize
    file '*trimming_report.txt' into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
    tg_length = "--length ${params.length}"
    c_r1 = params.clip_R1 > 0 ? "--clip_R1 ${params.clip_R1}" : ''
    tpc_r1 = params.three_prime_clip_R1 > 0 ? "--three_prime_clip_R1 ${params.three_prime_clip_R1}" : ''
    """
    trim_galore --small_rna $tg_length $c_r1 $tpc_r1 --gzip $reads --fastqc
    """
}


/*
 * STEP 2.1 - Insertsize
 */

process insertsize {
    tag "$reads"
    publishDir "${params.outdir}/trim_galore/insertsize", mode: 'copy'

    input:
    file reads from trimmed_reads_insertsize

    output:
    file '*.insertsize' into insertsize_results

    script:
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' <(zcat $reads) >${prefix}.insertsize
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

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    seqCenter = params.seqCenter ? "--sam-RG ID:${prefix} --sam-RG 'CN:${params.seqCenter}'" : ''
    """
    bowtie \\
        $index_base \\
        -q <(zcat $reads) \\
        -p 2 \\
        -t \\
        -k 1 \\
        -m 1 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        $seqCenter \\
        --un ${prefix}.mature_unmapped.fq \\
        -S \\
        | samtools view -bS - > ${prefix}.mature.bam

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

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - '.mature_unmapped.fq.gz'
    seqCenter = params.seqCenter ? "--sam-RG ID:${prefix} --sam-RG 'CN:${params.seqCenter}'" : ''
    """
    bowtie \\
        $index_base \\
        -p 2 \\
        -t \\
        -k 1 \\
        -m 1 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        -q <(zcat $reads) \\
        $seqCenter \\
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
    edgeR_miRBase.r $input_files
    """
}


/*
 * STEP 7.1 and 7.2 IF A GENOME SPECIFIED ONLY!
 */
if( params.gtf && params.bt_index) {

    /*
     * STEP 7.1 - Bowtie 1 against reference genome
     */
    process bowtie_ref {
        tag "$reads"
        publishDir "${params.outdir}/bowtie_ref", mode: 'copy'

        input:
        file reads from trimmed_reads_bowtie_ref
        file bt_indices

        output:
        file '*.bowtie.bam' into bowtie_bam, bowtie_bam_for_unmapped

        script:
        index_base = bt_indices[0].toString().tokenize(' ')[0].tokenize('.')[0]
        prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
        seqCenter = params.seqCenter ? "--sam-RG ID:${prefix} --sam-RG 'CN:${params.seqCenter}'" : ''
        """
        bowtie \\
            $index_base \\
            -q <(zcat $reads) \\
            -p 8 \\
            -t \\
            -k 10 \\
            -m 1 \\
            --best \\
            --strata \\
            -e 99999 \\
            --chunkmbs 2048 \\
            $seqCenter \\
            -S \\
            | samtools view -bS - > ${prefix}.bowtie.bam
        """
    }

    /*
     * STEP 7.2 - Statistics about unmapped reads against ref genome
     */

    process bowtie_unmapped {
        tag "${input_files[0].baseName}"
        publishDir "${params.outdir}/bowtie_ref/unmapped", mode: 'copy'

        input:
        file input_files from bowtie_bam_for_unmapped.toSortedList()

        output:
        file 'unmapped_refgenome.txt' into bowtie_unmapped

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
     * STEP 7.3 - NGI-Visualizations of Bowtie 1 alignment against host reference genome
     */
    process ngi_visualizations {
        tag "$bowtie_bam"
        publishDir "${params.outdir}/bowtie_ref/ngi_visualizations", mode: 'copy'

        input:
        file gtf from gtf
        file bowtie_bam

        output:
        file '*.{png,pdf}' into bowtie_ngi_visualizations

        script:
        // Note! ngi_visualizations needs to be installed!
        // See https://github.com/NationalGenomicsInfrastructure/ngi_visualizations
        """
        #!/usr/bin/env python
        from ngi_visualizations.biotypes import count_biotypes
        count_biotypes.main('$gtf','$bowtie_bam')
        """
    }

}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo "$params.version" > v_nfcore_smrnaseq.txt
    echo "$workflow.nextflow.version" > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version > v_trim_galore.txt
    bowtie --version > v_bowtie.txt
    samtools --version > v_samtools.txt
    fasta_formatter -h > v_fastx.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * STEP 8 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.toList()
    file ('trim_galore/*') from trimgalore_results.toList()
    file ('software_versions/*') from software_versions_yaml.toList()

    output:
    file '*multiqc_report.html' into multiqc_html
    file '*multiqc_data' into multiqc_data

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/smrnaseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/smrnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
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
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

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
          log.info "[nf-core/smrnaseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/smrnaseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    smrnaseqlogo = new File("$baseDir/assets/smrnaseq_logo.png").bytes.encodeBase64().toString()
    scilifelablogo = new File("$baseDir/assets/SciLifeLab_logo.png").bytes.encodeBase64().toString()
    ngilogo = new File("$baseDir/assets/NGI_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:smrnaseqlogo/, "data:image/png;base64,$smrnaseqlogo")
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

    log.info "[nf-core/smrnaseq] Pipeline Complete"

}
