# ![nf-core/smrnaseq](docs/images/smrnaseq_logo.png)

[![Build Status](https://travis-ci.com/nf-core/smrnaseq.svg?branch=master)](https://travis-ci.com/nf-core/smrnaseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/smrnaseq.svg)](https://hub.docker.com/r/nfcore/smrnaseq)

## Introduction
**nf-core/smrnaseq** is a bioinformatics best-practice analysis pipeline used for small RNA sequencing data.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

### Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Basic annotation
    * Alignment against miRBase mature miRNA ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
      * Post-alignment processing of miRBase mature miRNA counts ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
      * Analysis on miRBase mature miRNA counts ([`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html))
    * Alignment against miRBase hairpin for the unaligned reads in step 3 ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
       * Post-alignment processing of miRBase hairpin counts ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
       * Analysis on miRBase hairpin counts  ([`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html))

4. Alignment against host reference genome ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
    * Post-alignment processing of alignment against host reference genome ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5. Collapse reads ([`seqcsluter`](https://seqcluster.readthedocs.io/mirna_annotation.html#processing-of-reads))
6. Alignment agains miRBase hairpin ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
7. miRNA and isomiR annotation using the previous alignment ([`mirtop`](https://github.com/miRTop/mirtop))
8. miRNA quality control ([`mirtrace`](https://github.com/friedlanderlab/mirtrace))
9. Present QC for raw read, alignment, and expression results ([`MultiQC`](http://multiqc.info/))

### Documentation
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
nf-core/smrnaseq was originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels (@ewels), Chuan Wang (@chuan-wang) and Rickard Hammarén (@Hammarn). Updated by Lorena Pantano (@lpantano) from MIT.

## Citation
You can cite the `nf-core` pre-print as follows:  
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
