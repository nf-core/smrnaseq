# NGI - Small RNA-seq BP
Nextflow pipeline for small RNA sequencing best practice analysis at the NGI at SciLifeLab in Stockholm, Sweden

Written by Phil Ewels (@ewels), Rickard Hammarén (@Hammarn) and Chuan Wang (@chuan-wang)

# Under Development!
> This pipeline is currently being written and is
> not yet ready for production use.

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

## Running the pipeline
The typical command for running the pipeline is as follows:
```
nextflow run SciLifeLab/smNGI-RNAseq --reads '*_R{1,2}.fastq.gz'
```
or using a more manual approach (require you to clone the git repository)

```
nextflow path_to_NGI-smRNAseq/main.nf -c path_to_your_nextflow_config --reads '*_R{1,2}.fastq.gz' --genome 'GRCm38'
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
 --reads 'path/to/data/sample_*_{1,2}.fastq'
```

**NB: Must be enclosed in quotes!**

Note that the `{1,2}` parentheses are required to specify paired end data. Running `--reads '*.fastq'` will treat
all files as single end. The file path should be in quotation marks to prevent shell glob expansion.

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory.

### `--genome`
The reference genome to use of the analysis, needs to be one of the genome specified in the config file.
The human `GRCh37` genome is set as default.
```
--genome 'GRCm38'
```

### Supported genomes
```
Parameter     :       Species
AGPv3         :       Zea mays (Maize)
BDGP6         :       Drosophila melanogaster (Fruit fly)
CanFam3.1     :       Canis familiaris (Dog)
CHIMP2.1.4    :       Pan troglodytes (Chimpanze)
EquCab2       :       Equus caballus (Horse)
Galgal4       :       Gallus gallus (Chicken)
Gm01          :       Glycine max (Soybean)
GRCh37        :       Homo sapiens (Human)
GRCm38        :       Mus musculus (Mouse)
GRCz10        :       Danio rerio (Zebrafish)
IRGSP-1.0     :       Oryza sativa japonica (Rice)
Mmul_1        :       Macaca mulatta (Macaque)
Rnor_6.0      :       Rattus norvegicus (Rat)
Sbi1          :       Sorghum bicolor (Great millet)
Sscrofa10.2   :       Sus scrofa (Pig)
TAIR10        :       Arabidopsis thaliana (Thale cress)
UMD3.1        :       Bos taurus (Cow)
WBcel235      :       Caenorhabditis elegans (Nematode)
```

The `example_uppmax_config` file currently has the location of references for `GRCh37` (Human), `GRCm38` (Mouse)
and `sacCer2` (Yeast).

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes.
