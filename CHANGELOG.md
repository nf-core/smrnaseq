# nf-core/smrnaseq: Changelog

## v1.1dev - [date]

* Accept custom genome and remove non-canonical letters in the genome. Thanks to @sdjebali. Follow up from [[#63]](https://github.com/nf-core/smrnaseq/pull/63)
* Fix error when only one sample is in the input [[#31]](https://github.com/nf-core/smrnaseq/issues/31)
* Change `--reads` to `--input` for consistency with rest of nf-core
* Made `CamelCase` pipeline parameters `snake_case` and lower case
  * `clip_R1` -> `clip_r1`
  * `three_prime_clip_R1` -> `three_prime_clip_r1`
  * `saveReference` -> `save_reference`
  * `skipQC` -> `skip_qc`
  * `skipFastqc` -> `skip_fastqc`
  * `skipMultiqc` -> `skip_multiqc`
* Move all parameter documentation into new `nextflow_schema.json` file.
* Update with the latest `TEMPLATE` version for nf-core `1.12.1`
* Update conda environment with new packages and updates
* Added `--protocol custom` to allow custom adapter trimming modes [[#41]](https://github.com/nf-core/smrnaseq/issues/41)]
* Fix error when UMI is on the reads header: [[#35](https://github.com/nf-core/smrnaseq/issues/35)]
* Updated `params.mirtrace_species` to be properly initialised when using `--genome`, for all iGenomes species
* Made `params.mature` and `params.hairpin` default to miRBase FTP URL so that the file is automatically downloaded if not provided
* Allow `.fa` or `.fa.gz` files for mature and hairpin FASTA files.
* Add full-size benchmark / test dataset to run on AWS for each release (see `test_full.config`)
* Add `-f` flag to `gunzip` commands to deal with soft-links

### Packaged software updates

* `fastqc=0.11.8` -> `0.11.9`
* `trim-galore=0.6.3` -> `0.6.5`
* `bowtie=1.2.2` -> `1.2.3`
* `multiqc=1.7` -> `1.9`
* `mirtop=0.4.22` -> `0.4.23`
* `seqcluster=1.2.5` -> `1.2.7`
* `htseq=0.11.2` -> `0.11.3`
* `fastx_toolkit=0.0.14` -> `0.0.14`
* `seqkit=0.10.1` -> `0.12.0`
* `mirtrace=1.0.0` -> `1.0.1`
* Added `conda-forge::python=3.7.3`
* Added `conda-forge::markdown=3.1.1`
* Added `conda-forge::pymdown-extensions=6.0`
* Added `conda-forge::pygments=2.5.2`
* Removed `conda-forge::r-markdown=1.0`

## [v1.0.0](https://github.com/nf-core/smrnaseq/releases/tag/1.0.0) - 2019-09-19

### Added

* Figures to output documentation
* Samtools stats for genome alignments
* Seqkit and remove razers
* Mirtop and razers tools
* Adapt code and docs to [nf-core](http://nf-co.re/) template
* Update tools and Dockerfile/Singularity to match current template

### Packaged software updates

* openjdk 8.0.144 -> 11.0.1
* fastqc 0.11.7 -> 0.11.8
* trim-galore 0.5.0 -> 0.6.2
* bioconductor-edger 3.20.7 -> 3.26.0
* bioconductor-limma 3.34.9 -> 3.40.0
* conda-forge::r-data.table 1.11.4 -> 1.12.2
* conda-forge::r-gplots 3.0.1 -> 3.0.1.1
* conda-forge::r-r.methodss3 1.7.1 -> 1.7.1
* htseq 0.9.1 -> 0.11.2
* r-markdown 0.9
* Added mirtop 0.4.18a
* Removed razers3 3.5.3
* Added seqkit 0.10.1-1
* Added seqcluster 1.2.5
* conda-forge::r-base=3.5.1 -> 3.6.1
* conda-forge::r-statmod=1.4.30 -> 1.4.32
* conda-forge::r-markdown=0.9 -> 1.0
* trim-galore=0.6.2 -> 0.6.3
* mirtop=0.4.18a -> 0.4.22
* bioconductor-edger=3.26.0 -> 3.26.5
* bioconductor-limma=3.40.0 -> 3.40.2

## 2019-01-10

### Added

* "protocol" with pre-defined settings
* miRTrace in the pipeline.

### Software updates

* multiqc 1.6 -> 1.7.

## 2018-08-06

### Added

* Port original pipeline [SciLifeLab/NGI-smRNAseq](https://github.com/SciLifeLab/NGI-smRNAseq) to [nf-core/smrnaseq](https://github.com/nf-core/smrnaseq).
* Use Bowtie 1 instead of Bowtie 2 for the alignment to host reference genome.
* Option for sequencing centre in BAM file.

### Software updates

* trim-galore 0.4.5 -> 0.5.0
* samtools 1.8 -> 1.9
* multiqc 1.5 -> 1.6
