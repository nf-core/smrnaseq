# ![nf-core/smrnaseq](docs/images/nf-core-smrnaseq_logo_light.png#gh-light-mode-only) ![nf-core/smrnaseq](docs/images/nf-core-smrnaseq_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/smrnaseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/smrnaseq/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/smrnaseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/smrnaseq/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/smrnaseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3456879-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3456879)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/smrnaseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23smrnaseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/smrnaseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/smrnaseq** is a bioinformatics best-practice analysis pipeline for Small RNA-Seq.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/smrnaseq/results).

## Online videos

A short talk about the history, current status and functionality on offer in this pipeline was given by Lorena Pantano (@lpantano) on [9th November 2021](https://youtu.be/4YLQ2VwpCJE) as part of the nf-core/bytesize series.

You can find numerous talks on the nf-core events page from various topics including writing pipelines/modules in Nextflow DSL2, using nf-core tooling, running nf-core pipelines as well as more generic content like contributing to Github. Please check them out!

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
   1. Insert Size calculation
   2. Collapse reads ([`seqcluster`](https://seqcluster.readthedocs.io/mirna_annotation.html#processing-of-reads))
3. Contamination filtering ([`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
4. Alignment against miRBase mature miRNA ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
5. Alignment against miRBase hairpin
   1. Unaligned reads from step 3 ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
   2. Collapsed reads from step 2.2 ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
6. Post-alignment processing of miRBase hairpin
   1. Basic statistics from step 3 and step 4.1 ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
   2. Analysis on miRBase, or MirGeneDB hairpin counts ([`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html))
      - TMM normalization and a table of top expression hairpin
      - MDS plot clustering samples
      - Heatmap of sample similarities
   3. miRNA and isomiR annotation from step 4.1 ([`mirtop`](https://github.com/miRTop/mirtop))
7. Alignment against host reference genome ([`Bowtie1`](http://bowtie-bio.sourceforge.net/index.shtml))
   1. Post-alignment processing of alignment against host reference genome ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
8. Novel miRNAs and known miRNAs discovery ([`MiRDeep2`](https://www.mdc-berlin.de/content/mirdeep2-documentation))
   1. Mapping against reference genome with the mapper module
   2. Known and novel miRNA discovery with the mirdeep2 module
9. miRNA quality control ([`mirtrace`](https://github.com/friedlanderlab/mirtrace))
10. Present QC for raw read, alignment, and expression results ([`MultiQC`](http://multiqc.info/))

## Usage

:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1
Clone1_N1,s3://ngi-igenomes/test-data/smrnaseq/C1-N1-R1_S4_L001_R1_001.fastq.gz
Clone1_N3,s3://ngi-igenomes/test-data/smrnaseq/C1-N3-R1_S6_L001_R1_001.fastq.gz
Clone9_N1,s3://ngi-igenomes/test-data/smrnaseq/C9-N1-R1_S7_L001_R1_001.fastq.gz
Clone9_N2,s3://ngi-igenomes/test-data/smrnaseq/C9-N2-R1_S8_L001_R1_001.fastq.gz
Clone9_N3,s3://ngi-igenomes/test-data/smrnaseq/C9-N3-R1_S9_L001_R1_001.fastq.gz
Control_N1,s3://ngi-igenomes/test-data/smrnaseq/Ctl-N1-R1_S1_L001_R1_001.fastq.gz
Control_N2,s3://ngi-igenomes/test-data/smrnaseq/Ctl-N2-R1_S2_L001_R1_001.fastq.gz
Control_N3,s3://ngi-igenomes/test-data/smrnaseq/Ctl-N3-R1_S3_L001_R1_001.fastq.gz
```

Each row represents a fastq file (single-end).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/smrnaseq \
   -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --genome 'GRCh37' \
  --mirtrace_species 'hsa' \
  --protocol 'illumina' \
  --outdir <OUTDIR>
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/smrnaseq/usage) and the [parameter documentation](https://nf-co.re/smrnaseq/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/smrnaseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/smrnaseq/output).

## Credits

nf-core/smrnaseq was originally written by P. Ewels, C. Wang, R. Hammarén, L. Pantano, A. Peltzer.

We thank the following people for their extensive assistance in the development of this pipeline:

Lorena Pantano ([@lpantano](https://github.com/lpantano)) from MIT updated the pipeline to Nextflow DSL2.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#smrnaseq` channel](https://nfcore.slack.com/channels/smrnaseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/smrnaseq for your analysis, please cite it using the following doi: [10.5281/zenodo.3456879](https://zenodo.org/badge/latestdoi/140590861)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
