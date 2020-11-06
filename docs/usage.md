# nf-core/smrnaseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/smrnaseq/usage](https://nf-co.re/smrnaseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/smrnaseq --reads '*.fastq.gz' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/smrnaseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/smrnaseq releases page](https://github.com/nf-core/smrnaseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/smrnaseq`](https://hub.docker.com/r/nfcore/smrnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/smrnaseq`](https://hub.docker.com/r/nfcore/smrnaseq/)
* `podman`
  * A generic configuration profile to be used with [Podman](https://podman.io/)
  * Pulls software from Docker Hub: [`nfcore/smrnaseq`](https://hub.docker.com/r/nfcore/smrnaseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity or Podman.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
 --reads 'path/to/data/*.fastq.gz'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character

### `--protocol`

Protocol for constructing smRNA-seq libraries. Note that trimming parameters and 3' adapter sequence are pre-defined with a specified protocol.
Default: "illumina"

```bash
--protocol [one protocol listed in the table below]
```

| Protocol      | Library Prep Kit                        | Trimming Parameter                   | 3' Adapter Sequence   |
| :------------ | :-------------------------------------- | :----------------------------------- | :-------------------  |
| illumina      | Illumina TruSeq Small RNA               | clip_R1 = 0; three_prime_clip_R1 = 0 | TGGAATTCTCGGGTGCCAAGG |
| nextflex      | BIOO SCIENTIFIC  NEXTFLEX Small RNA-Seq | clip_R1 = 4; three_prime_clip_R1 = 4 | TGGAATTCTCGGGTGCCAAGG |
| qiaseq        | QIAGEN QIAseq miRNA                     | clip_R1 = 0; three_prime_clip_R1 = 0 | AACTGTAGGCACCATCAAT   |
| cats          | Diagenode CATS Small RNA-seq            | clip_R1 = 3; three_prime_clip_R1 = 0 | AAAAAAAAAAA + GATCGGAAGAGCACACGTCTG (only polyA is used for trimming) |

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
      mature  = '<path to the genome fasta file>' //mature.fa"
      hairpin = '<path to the genome fasta file>' //hairpin.fa"
      bowtie  = '<path to the genome fasta file>' // index
      gtf     = '<path to the gtf file>'
      mirtrace_species = "sps" // species according mirbase

    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--saveReference`

Supply this parameter to save any generated reference genome files to your results folder. These can then be used for future pipeline runs, reducing processing times.

### `--fasta`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

### `--igenomes_ignore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

### `--mature`

If you prefer, you can specify the full path to the FASTA file of mature miRNAs when you run the pipeline:

```bash
--mature [path to the FASTA file of mature miRNAs]
```

### `--hairpin`

If you prefer, you can specify the full path to the FASTA file of miRNA precursors when you run the pipeline:

```bash
--hairpin [path to the FASTA file of miRNA precursors]
```

### `--bt_index`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--bt_index [path to Bowtie 1 index]
```

## Trimming options

### `--min_length [int]`

Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18

### `--clip_R1 [int]`

Instructs Trim Galore to remove bp from the 5' end of read 1

### `--three_prime_clip_R1 [int]`

Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed

### `--three_prime_adapter [sequence]`

Instructs Trim Galore to remove 3' adapters which are typically used in smRNA-seq library preparation

## Skipping QC steps

### `--skipQC`

Skip all QC steps aside from MultiQC

### `--skipFastqc`

Skip FastQC

### `--skipMultiqc`

Skip MultiQC

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `--seq_center`

Text about sequencing center which will be added in the header of output bam files.

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, we run on UPPMAX but don't want to use the MultiQC
environment module as is the default. So we specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

## Stand-alone scripts

The `bin` directory contains some scripts used by the pipeline which may also be run manually:

* `edgeR_miRBase.r`
  * R script using for processing reads counts of mature miRNAs and miRNA precursors (hairpins).

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
