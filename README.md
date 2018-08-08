# ![nf-core/smrnaseq](docs/images/smrnaseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/smrnaseq.svg?branch=master)](https://travis-ci.org/nf-core/smrnaseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.2-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/Lobby)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/smrnaseq.svg)](https://hub.docker.com/r/nfcore/smrnaseq/)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1250)


----

# UNDER DEVELOPMENT!
This pipeline has recently been moved to nf-core and is still under heavy development. It does not yet meet all of the requirements for nf-core pipelines.

Use with caution!

----

**nf-core/smrnaseq** is a bioinformatics best-practice analysis pipeline used for small RNA sequencing data. It is developed at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This pipeline is primarily used with a SLURM cluster on the Swedish [UPPMAX systems](https://www.uppmax.uu.se). However, the pipeline should be able to run on any system that Nextflow supports. We have done some limited testing using Docker and AWS, and the pipeline comes with some configuration for these systems. See the [installation docs](docs/installation.md) for more information.

## Installation
### NextFlow installation
See https://www.nextflow.io for instructions on how to install and configure
Nextflow.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if
`nf-core/smrnaseq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:
```
git clone https://github.com/nf-core/smrnaseq.git
nextflow run nf-core/smrnaseq/main.nf
```

## Configuration
For running this pipeline on the the Swedish UPPMAX cluster, command line flag `--project <project_ID>`.

To avoid having to specify this every time you run Nextflow, you can add it to your
personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```nextflow
params.project = 'project_ID'
```

The pipeline will exit with an error message if you try to run it pipeline with the default
UPPMAX config profile and don't set project.


### Running on other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up
your own config file so that the script knows where to find your reference files and how your
cluster works.

Copy the contents of [`conf/uppmax.config`](conf/uppmax.config) to your own config file somewhere
and then reference it with `-c` when running the pipeline.

If you think that there are other people using the pipeline who would benefit from your configuration
(eg. other common cluster setups), please let us know. It should be easy to create a new config file
in `conf` and reference this as a named profile in [`nextflow.config`](nextflow.config). Then these
configuration options can be used by specifying `-profile <name>` when running the pipeline.


## Running the pipeline
The typical command for running the pipeline is as follows:

```
nextflow run nf-core/smrnaseq --reads '*.fastq.gz'
```

**NOTE! Paired-end data is NOT supported by this pipeline!**
For paired-end data, use Read 1 only. For instance:

```
nextflow run nf-core/smrnaseq --reads '*.R1.fastq.gz'
```

Note that the pipeline will create files in your working directory:
```bash
work            # Directory containing the nextflow working files
results         # Finished results for each sample, one directory per pipeline step
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Mandatory parameters
### `--reads`
Location of the input FastQ files:
```
 --reads 'path/to/data/*.fastq.gz'
```
**NOTE! Must be enclosed in quotes!**
If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory.

### `--genome`
The reference genome to use of the analysis, needs to be one of the genome specified in the config file.
The human `GRCh37` genome is used by default.
```
--genome 'GRCh37'
```

### Supported genomes   

| Parameter     |       Latin Name                 |      Common Name   |
| :------------ |:-------------------------------- |:------------------ |
| AGPv3         |       *Zea mays*                 |       Maize        |
| BDGP6         |       *Drosophila melanogaster*  |       Fruit fly    |
| CanFam3.1     |       *Canis familiaris*         |       Dog          |
| CHIMP2.1.4    |       *Pan troglodytes*          |       Chimpanze    |
| EquCab2       |       *Equus caballus*           |       Horse        |
| Galgal4       |       *Gallus gallus*            |       Chicken      |
| Gm01          |       *Glycine max*              |       Soybean      |
| GRCh37        |       *Homo sapiens*             |       Human        |
| GRCm38        |       *Mus musculus*             |       Mouse        |
| GRCz10        |       *Danio rerio*              |       Zebrafish    |
| IRGSP-1.0     |       *Oryza sativa japonica*    |       Rice         |
| Mmul_1        |       *Macaca mulatta*           |       Macaque      |
| Rnor_6.0      |       *Rattus norvegicus*        |       Rat          |
| Sbi1          |       *Sorghum bicolor*          |       Great millet |
| Sscrofa10.2   |       *Sus scrofa*               |       Pig          |
| TAIR10        |       *Arabidopsis thaliana*     |       Thale cress  |
| UMD3.1        |       *Bos taurus*               |       Cow          |
| WBcel235      |       *Caenorhabditis elegans*   |       Nematode     |


## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `--seqCenter`
Text about sequencing center which will be added in the header of output bam files. Note that no blank is allowed!

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes. **NOTE! One hyphen only (core Nextflow parameter).**

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, we run on UPPMAX but don't want to use the MultiQC
environment module as is the default. So we specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

### `--mature`
If you prefer, you can specify the full path to the FASTA file of mature miRNAs when you run the pipeline:
```
--mature [path to the FASTA file of mature miRNAs]
```

### `--hairpin`
If you prefer, you can specify the full path to the FASTA file of miRNA precursors when you run the pipeline:
```
--hairpin [path to the FASTA file of miRNA precursors]
```

### `--bt_index`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:
```
--bt_index [path to Bowtie 1 index]
```

### `--rlocation`
Some steps in the pipeline run R with required modules. By default, the pipeline will install
these modules to `~/R/nxtflow_libs/` if not present. You can specify what path to use with this
command line flag.

### Trimming options
`--length [int]`: Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18
`--clip_R1 [int]`: Instructs Trim Galore to remove bp from the 5' end of read 1
`--three_prime_clip_R1 [int]`: Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder. These can then be used for future pipeline runs, reducing processing times.

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used instead of the config file specific to the pipeline.

### `--clusterOptions`
Submit arbitrary SLURM options (UPPMAX profile only). For instance, you could use `--clusterOptions '-p devcore'`
to run on the development node (though won't work with default process time requests).

## Stand-alone scripts
The `bin` directory contains some scripts used by the pipeline which may also be run manually:

* `edgeR_miRBase.r`
  * R script using for processing reads counts of mature miRNAs and miRNA precursors (hairpins).


## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

Written by Phil Ewels (@ewels), Chuan Wang (@chuan-wang) and Rickard Hammarén (@Hammarn)

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a></p>
