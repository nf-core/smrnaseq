# nf-core/smrnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2.1.0](https://github.com/nf-core/smrnaseq/releases/tag/2.1.0) - 2022-10-25 Maroon Tin Dalmatian

### Enhancements & fixes

- [[#12](https://github.com/nf-core/smrnaseq/issues/12)] - Enabled the use of `MirGeneDB` as an alternative database insted of `miRBase`
- [[#113](https://github.com/nf-core/smrnaseq/issues/113)] - Added a optional contamination filtering step, including MultiQC plot
- [[#137](https://github.com/nf-core/smrnaseq/issues/137)] - Fixed issue with mirTop and MultiQC by upgrading to MultiQC V1.13dev
- [[#159](https://github.com/nf-core/smrnaseq/issues/159)] - Index files were not collected when `bowtie_index` was used and thus mapping was failing
- [[#161](https://github.com/nf-core/smrnaseq/issues/161)] - Trimmed output was not as documented and not correctly published
- [[#168](https://github.com/nf-core/smrnaseq/issues/168)] - Removed `mirtrace_protocol` as the parameter was redundant and `params.protocol` is entirely sufficient
- Updated pipeline template to [nf-core/tools 2.6.0](https://github.com/nf-core/tools/releases/tag/2.6.0)
- [[#188](https://github.com/nf-core/smrnaseq/pull/188)] - Dropped TrimGalore in favor of fastp QC and adapter trimming, improved handling of adapters and trimming parameters
- [[#194](https://github.com/nf-core/smrnaseq/issues/194)] - Added default adapters file for FastP improved miRNA adapter trimming

### Parameters

| Old parameter | New parameter            |
| ------------- | ------------------------ |
|               | `--mirgenedb`            |
|               | `--mirgenedb_species`    |
|               | `--mirgenedb_gff`        |
|               | `--mirgenedb_mature`     |
|               | `--mirgenedb_hairpin`    |
|               | `--contamination_filter` |
|               | `--rrna`                 |
|               | `--trna`                 |
|               | `--cdna`                 |
|               | `--ncrna`                |
|               | `--pirna`                |
|               | `--other_contamination`  |

## [v2.0.0](https://github.com/nf-core/smrnaseq/releases/tag/2.0.0) - 2022-05-31 Aqua Zinc Chihuahua

### Major enhancements

- Updated pipeline template to [nf-core/tools 2.4.1](https://github.com/nf-core/tools/releases/tag/2.4.1)
- [[#137](https://github.com/nf-core/smrnaseq/issues/137)] - Update mirtop container to version `0.4.25` to fix multiqc error
- Port pipeline to the updated Nextflow DSL2 syntax adopted on nf-core/modules
- Bump minimum Nextflow version from `20.04.0` -> `21.10.3`
- Point to the proper test data branch
- Software version(s) will now be reported for every module imported during a given pipeline execution
- Updated `nextflow_schema.json` should now display correctly on Nextflow Tower
- Added mirtop logs to multiqc
- Allow a gene to be associated to a non null number of reads in all samples (in `edgeR_miRBase.r` script)

### Other enhancements & fixes

- [#134](https://github.com/nf-core/smrnaseq/issues/134) - Fixed colSum of zero issues for edgeR_miRBase.R script
- [#55](https://github.com/lpantano/seqcluster/pull/55) - update seqcluster to fix UMI-detecting bug

### Parameters

| Old parameter        | New parameter    |
| -------------------- | ---------------- |
| `--conda`            | `--enable_conda` |
| `--clusterOptions`   |                  |
| `--publish_dir_mode` |                  |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency           | Old version | New version |
| -------------------- | ----------- | ----------- |
| `bioconductor-edger` | 3.26.5      | 3.36.0      |
| `bioconductor-limma` | 3.40.2      | 3.50.0      |
| `mirdeep2`           | 2.0.1.2     | 2.0.1.3     |
| `mirtop`             | 0.4.23      | 0.4.25      |
| `multiqc`            | 1.10.1      | 1.12.0      |
| `python`             | 3.7.3       | 3.8.3       |
| `r-base`             | 3.6.3       | 4.0.3       |
| `r-data.table`       | 1.12.4      | 1.14.2      |
| `r-gplots`           | 3.0.1.1     | 3.1.1       |
| `r-statmod`          | 1.4.32      | 1.4.36      |
| `samtools`           | 1.12        | 1.13        |
| `seqcluster`         | 1.2.7       | 1.2.9       |
| `seqkit`             | 0.16.0      | 2.0.0       |
| `trim-galore`        | 0.6.6       | 0.6.7       |
| `bioconvert`         | -           | 0.4.3       |
| `htseq`              | -           | -           |
| `markdown`           | -           | -           |
| `pymdown-extensions` | -           | -           |
| `pygments`           | -           | -           |
| `r-r.methodss3`      | -           | -           |
| `bowtie2`            | -           | 2.4.5       |
| `blat`               | -           | 36          |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [v1.1.0](https://github.com/nf-core/smrnaseq/releases/tag/1.1.0) - 2021-06-15

### Major changes

**:warning: Breaking changes!**

This release contains several major (potentially breaking) changes:

- The main input parameter has been changed from `--reads` to `--input` to standardize the pipeline with other nf-core pipelines
- All parameter documentation has moved into a new `nextflow_schema.json` file
- A `lib` folder with groovy helper classes has been added to the pipeline. This includes validation of input parameters using the schema defined in the `nextflow_schema.json` file

### General improvements

- remove spaces in genome headers and replace special nt by N in hairpin file for mirdeep2 to work. [[#69]](https://github.com/nf-core/smrnaseq/pull/79)
- Accept custom genome and remove non-canonical letters in the genome. Thanks to @sdjebali. Follow up from [[#63]](https://github.com/nf-core/smrnaseq/pull/63)
- Fix error when only one sample is in the input [[#31]](https://github.com/nf-core/smrnaseq/issues/31)
- Made `CamelCase` pipeline parameters `snake_case` and lower case
  - `clip_R1` -> `clip_r1`
  - `three_prime_clip_R1` -> `three_prime_clip_r1`
  - `saveReference` -> `save_reference`
  - `skipQC` -> `skip_qc`
  - `skipFastqc` -> `skip_fastqc`
  - `skipMultiqc` -> `skip_multiqc`
- Update with the latest `TEMPLATE` version for nf-core `1.12.1`
- Update conda environment with new packages and updates
- Added `--protocol custom` to allow custom adapter trimming modes [[#41]](https://github.com/nf-core/smrnaseq/issues/41)]
- Fix error when UMI is on the reads header: [[#35](https://github.com/nf-core/smrnaseq/issues/35)]
- Updated `params.mirtrace_species` to be properly initialised when using `--genome`, for all iGenomes species
- Made `params.mature` and `params.hairpin` default to miRBase FTP URL so that the file is automatically downloaded if not provided
- Allow `.fa` or `.fa.gz` files for mature and hairpin FASTA files.
- Add full-size benchmark / test dataset to run on AWS for each release (see `test_full.config`)
- Add `-f` flag to `gunzip` commands to deal with soft-links
- Add `--mirtrace_protocol` parameter to allow for more flexible selection of this parameter
- Added `--trim_galore_max_length` option [[#77](https://github.com/nf-core/smrnaseq/issues/77)]
- Solved memory usage issue for mirtrace in the main process and in the `get_software_versions` process [[#68](https://github.com/nf-core/smrnaseq/issues/68)]
- Removed logging of `single_end` parameter and added missing parameters to schema and config files
- Added "custom" as option for `--protocol` in the `nextflow_schema.json`

### Packaged software updates

- `fastqc=0.11.8` -> `0.11.9`
- `trim-galore=0.6.3` -> `0.6.5`
- `bowtie=1.2.2` -> `1.2.3`
- `multiqc=1.7` -> `1.9`
- `mirtop=0.4.22` -> `0.4.23`
- `seqcluster=1.2.5` -> `1.2.7`
- `htseq=0.11.2` -> `0.11.3`
- `fastx_toolkit=0.0.14` -> `0.0.14`
- `seqkit=0.10.1` -> `0.12.0`
- `mirtrace=1.0.0` -> `1.0.1`
- Added `conda-forge::python=3.7.3`
- Added `conda-forge::markdown=3.1.1`
- Added `conda-forge::pymdown-extensions=6.0`
- Added `conda-forge::pygments=2.5.2`
- Removed `conda-forge::r-markdown=1.0`

## [v1.0.0](https://github.com/nf-core/smrnaseq/releases/tag/1.0.0) - 2019-09-19

### Added

- Figures to output documentation
- Samtools stats for genome alignments
- Seqkit and remove razers
- Mirtop and razers tools
- Adapt code and docs to [nf-core](http://nf-co.re/) template
- Update tools and Dockerfile/Singularity to match current template

### Packaged software updates

- openjdk 8.0.144 -> 11.0.1
- fastqc 0.11.7 -> 0.11.8
- trim-galore 0.5.0 -> 0.6.2
- bioconductor-edger 3.20.7 -> 3.26.0
- bioconductor-limma 3.34.9 -> 3.40.0
- conda-forge::r-data.table 1.11.4 -> 1.12.2
- conda-forge::r-gplots 3.0.1 -> 3.0.1.1
- conda-forge::r-r.methodss3 1.7.1 -> 1.7.1
- htseq 0.9.1 -> 0.11.2
- r-markdown 0.9
- Added mirtop 0.4.18a
- Removed razers3 3.5.3
- Added seqkit 0.10.1-1
- Added seqcluster 1.2.5
- conda-forge::r-base=3.5.1 -> 3.6.1
- conda-forge::r-statmod=1.4.30 -> 1.4.32
- conda-forge::r-markdown=0.9 -> 1.0
- trim-galore=0.6.2 -> 0.6.3
- mirtop=0.4.18a -> 0.4.22
- bioconductor-edger=3.26.0 -> 3.26.5
- bioconductor-limma=3.40.0 -> 3.40.2

## 2019-01-10

### Added

- "protocol" with pre-defined settings
- miRTrace in the pipeline.

### Software updates

- multiqc 1.6 -> 1.7.

## 2018-08-06

### Added

- Port original pipeline [SciLifeLab/NGI-smRNAseq](https://github.com/SciLifeLab/NGI-smRNAseq) to [nf-core/smrnaseq](https://github.com/nf-core/smrnaseq).
- Use Bowtie 1 instead of Bowtie 2 for the alignment to host reference genome.
- Option for sequencing centre in BAM file.

### Software updates

- trim-galore 0.4.5 -> 0.5.0
- samtools 1.8 -> 1.9
- multiqc 1.5 -> 1.6
