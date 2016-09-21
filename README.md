# NGI - Small RNA-seq BP
Nextflow pipeline for small RNA sequencing best practice analysis at the NGI at SciLifeLab in Stockholm, Sweden

Written by Phil Ewels (@ewels), Chuan Wang (@chuan-wang) and Rickard Hammarén (@Hammarn)


## Installation
### NextFlow installation
To use this pipeline, you need to have a working version of NextFlow installed. You can find more
information about this pipeline tool at [nextflow.io](http://www.nextflow.io/). The typical installation
of NextFlow looks like this:

```
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```
Note that if you're running on the Swedish UPPMAX cluster (Milou) you can load NextFlow as an
environment module:
```
module load nextflow
```

### NextFlow configuration
Next, you need to set up a config file so that NextFlow knows how to run and where to find reference
indexes. You can find an example configuration file for UPPMAX (milou) with this repository:
[`example_uppmax_config`](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/example_uppmax_config).

Copy this file to `~/.nextflow/config` and edit the line `'-A <PROJECT>'` to contain your own UPPMAX project
identifier instead.

It is entirely possible to run this pipeline on other clusters - just note that you may need to customise
the `process` environment (eg. if you're using a cluster system other than SLURM) and the paths to reference
files.

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

### Installation of the NGI plugin for the'MultiQC' module
```
pip install git+https://github.com/ewels/MultiQC_NGI.git
```

### IMPORTANT!
Currently there two steps that require local script files to run:  
1) Trim Galore!: We are using the latest release v0.4.2, in which a bug for small RNA-seq adapter trimming is fixed. This release is not installed in Uppmax yet.  
2) Note that for ngi_visualizations, python packages HTSeq and pysam are required.  

## Running the pipeline
The typical command for running the pipeline is as follows:

```
nextflow run SciLifeLab/smNGI-RNAseq --reads '*.fastq.gz'
```
or using a more manual approach (require you to clone the git repository)

```
nextflow path_to_NGI-smRNAseq/main.nf -c path_to_your_nextflow_config --reads '*.fastq.gz' --genome 'GRCh37'
```

**NOTE! Paired-end data is NOT supported by this pipeline!**
For paired-end data, use Read 1 only. For instance:

```
nextflow run SciLifeLab/smNGI-RNAseq --reads '*.R1.fastq.gz'
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
 --reads 'path/to/data/*.fastq'
```

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory.

### `--genome`
The reference genome to use of the analysis, needs to be one of the genome specified in the config file.
The human `GRCh37` genome is set as default.
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


### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes.
