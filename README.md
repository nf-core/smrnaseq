# ![nf-core/smrnaseq](docs/images/smrnaseq_logo.png)

[![Build Status](https://travis-ci.com/nf-core/smrnaseq.svg?branch=master)](https://travis-ci.com/nf-core/smrnaseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/smrnaseq.svg)](https://hub.docker.com/r/nfcore/smrnaseq)


## Pipeline summary:

- 1:   Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- 2:   Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
- 3.1: Alignment against miRBase mature miRNA ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
- 3.2: Post-alignment processing of miRBase mature miRNA counts ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
- 3.3: Analysis on miRBase mature miRNA counts ([`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html))
       - TMM normalization and a table of top expression mature miRNA
       - MDS plot clustering samples
       - Heatmap of sample similarities
- 4.1: Alignment against miRBase hairpin for the unaligned reads in step 3 ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
- 4.2: Post-alignment processing of miRBase hairpin counts ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
- 4.3: Analysis on miRBase hairpin counts  ([`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html))
       - TMM normalization and a table of top expression hairpin
       - MDS plot clustering samples
       - Heatmap of sample similarities
- 5.1: Alignment against host reference genome ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
- 5.2: Post-alignment processing of alignment against host reference genome ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
- 6:   Visualization of alignment statistics ([`NGI-Visualization`](https://github.com/NationalGenomicsInfrastructure/ngi_visualizations))
- 7:   miRNA quality control ([`mirtrace`](https://github.com/friedlanderlab/mirtrace))
- 8:  Present QC for raw read, alignment, and expression results ([`MultiQC`](http://multiqc.info/))

## Introduction

**nf-core/smrnaseq** is a bioinformatics best-practice analysis pipeline used for small RNA sequencing data.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Documentation
The nf-core/smrnaseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits
nf-core/smrnaseq was originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels (@ewels), Chuan Wang (@chuan-wang) and Rickard Hammarén (@Hammarn).

It is been updated by Lorena Pantano (@lpantano).
