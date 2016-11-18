# NGI - Small RNA-seq BP
Nextflow pipeline for small RNA sequencing best practice analysis at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

# UNDER DEVELOPMENT
This pipeline is currently under active development and has very little in the way of testing. You may well have problems running it. Any help debugging is very welcome! Please fork, make changes and submit a pull request.

## Installation
### NextFlow installation
See https://github.com/SciLifeLab/NGI-NextflowDocs for instructions on how to install and configure
Nextflow.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if
`SciLifeLab/NGI-smRNAseq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:
```
git clone https://github.com/SciLifeLab/NGI-smRNAseq.git
nextflow run NGI-smRNAseq/main.nf
```

### Installation of the 'ngi_visualizations' module
This module needs to be installed locally in order to visualize the statistics from Bowtie2 alignment.
```
pip install -U git+https://github.com/NationalGenomicsInfrastructure/ngi_visualizations.git
```
Note that for ngi_visualizations, python packages HTSeq and pysam are required.

### Installation of the NGI plugin for the'MultiQC' module
```
pip install git+https://github.com/ewels/MultiQC_NGI.git
```

## Configuration
By default, the pipeline is configured to run on the Swedish UPPMAX cluster (milou / irma).

You will need to specify your UPPMAX project ID when running a pipeline. To do this, use
the command line flag `--project <project_ID>`.

To avoid having to specify this every time you run Nextflow, you can add it to your
personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
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
nextflow run SciLifeLab/NGI-smRNAseq --reads '*.fastq.gz'
```

**NOTE! Paired-end data is NOT supported by this pipeline!**
For paired-end data, use Read 1 only. For instance:

```
nextflow run SciLifeLab/NGI-smRNAseq --reads '*.R1.fastq.gz'
```

Note that the pipeline will create files in your working directory:
```bash
work            # Directory containing the nextflow working files
results         # Finished results for each sample, one directory per pipeline step
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

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

**NOTE! With the option --genome 'ALL', the entire dataset of mature miRNAs and hairpins in miRBase will be used as reference regardless of species. Meanwhile the alignment against host reference genome will be skipped.**

### `--bt2index`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:
```
--bt2index [path to Bowtie2 index]
```

### `--outdir`
The output directory where the results will be saved.

### `--rlocation`
Some steps in the pipeline run R with required modules. By default, the pipeline will install
these modules to `~/R/nxtflow_libs/` if not present. You can specify what path to use with this
command line flag.

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes. **NOTE! One hyphen only (core Nextflow parameter).**


## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

Written by Phil Ewels (@ewels), Chuan Wang (@chuan-wang) and Rickard Hammarén (@Hammarn)

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a></p>
