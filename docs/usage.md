# nf-core/smrnaseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/smrnaseq/usage](https://nf-co.re/smrnaseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

### `protocol`

This option indicates the experimental protocol used for the sample preparation. Currently supporting:

- 'illumina': adapter (`TGGAATTCTCGGGTGCCAAGG`)
- 'nextflex': adapter (`TGGAATTCTCGGGTGCCAAGG), clip_r1 (`4`), three_prime_clip_r1 (`4`)
- 'qiaseq': adapter (`AACTGTAGGCACCATCAAT`)
- 'cats': adapter (`GATCGGAAGAGCACACGTCTG), clip_r1(`3)
- 'custom' (where the user can indicate the `three_prime_adapter`, `clip_r1` and `three_prime_clip_r1` manually)

:warning: At least the `custom` protocol has to be specified, otherwise the pipeline won't run. In case you specify the `custom` protocol, ensure that the parameters above are set accordingly or the defaults will be applied. If you want to auto-detect the adapters using `fastp`, please set `--three_prime_adapter` to `""`.

### `mirtrace_species` or `mirgenedb_species`

It should point to the 3-letter species name used by [miRBase](https://www.mirbase.org/help/genome_summary.shtml) or [MirGeneDB](https://www.mirgenedb.org/browse). Note the difference in case for the two databases.

### miRNA related files

Different parameters can be set for the two supported databases. By default `miRBase` will be used with the parameters below.

- `mirna_gtf`: If not supplied by the user, then `mirna_gtf` will point to the latest GFF3 file in miRbase: `https://mirbase.org/download/CURRENT/genomes/${params.mirtrace_species}.gff3`
- `mature`: points to the FASTA file of mature miRNA sequences. `https://mirbase.org/download/mature.fa`
- `hairpin`: points to the FASTA file of precursor miRNA sequences. `https://mirbase.org/download/hairpin.fa`

If MirGeneDB should be used instead it needs to be specified using `--mirgenedb` and use the parameters below.

- `mirgenedb_gff`: The data can not be downloaded automatically (URLs are created with short term tokens in it), thus the user needs to supply the gff file for either his species, or all species downloaded from `https://mirgenedb.org/download`. The total set will automatically be subsetted to the species specified with `--mirgenedb_species`.
- `mirgenedb_mature`: points to the FASTA file of mature miRNA sequences. Download from `https://mirgenedb.org/download`.
- `mirgenedb_hairpin`: points to the FASTA file of precursor miRNA sequences. Download from `https://mirgenedb.org/download`. Note that MirGeneDB does not have a dedicated `hairpin` file, but the `Precursor sequences` are to be used.

### Genome

- `fasta`: the reference genome FASTA file
- `bt_indices`: points to the folder containing the `bowtie2` indices for the genome reference specified by `fasta`. **Note:** if the FASTA file in `fasta` is not the same file used to generate the `bowtie2` indices, then the pipeline will fail.

### Contamination filtering

This step has, until now, only been tested for human data. Unexpected behaviour can occur when using it with a different species.

Contamination filtering of the sequencing reads is optional and can be invoked using the `filter_contamination` parameter. FASTA files with

- `rrna`: Used to supply a FASTA file containing rRNA contamination sequence.
- `trna`: Used to supply a FASTA file containing tRNA contamination sequence. e.g. `http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa`
- `cdna`: Used to supply a FASTA file containing cDNA contamination sequence. e.g. `ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz` The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `ncrna`: Used to supply a FASTA file containing ncRNA contamination sequence. e.g. `ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz` The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `pirna`: Used to supply a FASTA file containing piRNA contamination sequence. e.g. The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `other_contamination`: Used to supply an additional filtering set. The FASTA file is first compared to the available miRNA sequences and overlaps are removed.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 2 columns ("sample" and "fastq_1"), and a header row as shown in the examples below.

If a second fastq file is provided using another column, the extra data are ignored by this pipeline. The smRNA species should be sufficiently contained in the first read, and so the second read is superfluous data in this smRNA context.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers should match between runs of resequenced samples. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire. However, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/smrnaseq --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/smrnaseq -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/smrnaseq
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/smrnaseq releases page](https://github.com/nf-core/smrnaseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Stand-alone scripts

The `bin` directory contains some scripts used by the pipeline which may also be run manually:

- `edgeR_miRBase.r`: R script using for processing reads counts of mature miRNAs and miRNA precursors (hairpins).
  This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config once. However, if multiple users within an organisation intend to run nf-core pipelines regularly under the same settings, then consider uploading your custom config file to the `nf-core/configs` git repository. Before uploading, ensure that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amended [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
