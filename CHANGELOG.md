# nf-core/smrnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.4.0dev - 2024-XX-XX - X

- [[#349]](https://github.com/nf-core/smrnaseq/pull/349) - Fix [MIRTOP_QUANT conda issue](https://github.com/nf-core/smrnaseq/issues/347) - change conda-base to conda-forge channel.
- [[#350]](https://github.com/nf-core/smrnaseq/pull/350) - Fix [MIRTOP_QUANT conda issue](https://github.com/nf-core/smrnaseq/issues/347) - set python version to 3.7 to fix pysam issue.
- [[#361]](https://github.com/nf-core/smrnaseq/pull/361) - Fix [[#332]](https://github.com/nf-core/smrnaseq/issues/332) - Fix documentation to use only single-end.
- [[#364]](https://github.com/nf-core/smrnaseq/pull/364) - Fix [Protocol inheritance issue](https://github.com/nf-core/smrnaseq/issues/351) - fixing protocol inheritance from subworkflow with move to config profile(s) for different protocols.
- [[#372]](https://github.com/nf-core/smrnaseq/pull/372) - Fix [Plain test profile](https://github.com/nf-core/smrnaseq/issues/371) - Updated default protocol value to "custom".
- [[#374]](https://github.com/nf-core/smrnaseq/pull/374) - Fix [default tests](https://github.com/nf-core/smrnaseq/issues/375) so that they do not require additional profiles in CI. Change GitHub CI fail-fast strategy to false.
- [[#375]](https://github.com/nf-core/smrnaseq/pull/375) - Test [technical repeats](https://github.com/nf-core/smrnaseq/issues/212) - Test merging of technical repeats.
- [[#377]](https://github.com/nf-core/smrnaseq/pull/377) - Fix [Linting](https://github.com/nf-core/smrnaseq/issues/369) - Fixed linting warnings and updated modules & subworkflows.
- [[#378]](https://github.com/nf-core/smrnaseq/pull/378) - Fix [`--mirtrace_species` bug](<(https://github.com/nf-core/smrnaseq/issues/348)>) - Make `MIRTRACE` process conditional. Add mirgenedb test.
- [[#380]](https://github.com/nf-core/smrnaseq/pull/380) - Fix [edgeR_mirBase.R](https://github.com/nf-core/smrnaseq/issues/187) - Fix checking number of samples which causes error in plotMDS. Add nf-tests for local modules using custom R scripts.
- [[#381]](https://github.com/nf-core/smrnaseq/pull/381) - Update [Convert tests to nf-tests](https://github.com/nf-core/smrnaseq/issues/379) - CI tests to nf-tests.
- [[#382]](https://github.com/nf-core/smrnaseq/pull/382) - Add [collapse_mirtop.R](https://github.com/nf-core/smrnaseq/issues/174) - Add nf-tests for local modules using custom R scripts.
- [[#383]](https://github.com/nf-core/smrnaseq/pull/383) - Fix [parameter `--skip_fastp` throws an error](https://github.com/nf-core/smrnaseq/issues/263) - Fix parameter --skip_fastp.
- [[#384]](https://github.com/nf-core/smrnaseq/pull/384) - Fix [filter status bug fix](https://github.com/nf-core/smrnaseq/issues/360) - Fix filter stats module and add filter contaminants test profile.
- [[#386]](https://github.com/nf-core/smrnaseq/pull/386) - Fix [Nextflex trimming support](https://github.com/nf-core/smrnaseq/issues/365) - Fix Nextflex trimming support.
- [[#387]](https://github.com/nf-core/smrnaseq/pull/387) - Add [contaminant filter failure because the Docker image for BLAT cannot be pulled](https://github.com/nf-core/smrnaseq/issues/354) - Add nf-test to local module `blat_mirna` and fixes . Adds a small test profile to test contaminant filter results.
- [[#388]](https://github.com/nf-core/smrnaseq/pull/388) - Fix [igenomes fix](https://github.com/nf-core/smrnaseq/issues/360) - Fix workflow scripts so that they can use igenome parameters.
- [[#391]](https://github.com/nf-core/smrnaseq/pull/391) - Fix [error because of large chromosomes](https://github.com/nf-core/smrnaseq/issues/132) - Change `.bai` index for `.csi` index in `samtools_index` to fix .
- [[#392]](https://github.com/nf-core/smrnaseq/pull/392) - Update [Reduce tests](https://github.com/nf-core/smrnaseq/issues/389) - Combine and optimize tests, and reduce samplesheets sizes.
- [[#397]](https://github.com/nf-core/smrnaseq/pull/397) - Fix [contaminant filter failure because of the Docker image for BLAT](https://github.com/nf-core/smrnaseq/issues/354) - Improvements to contaminant filter subworkflow and replacement for nf-core modules.
- [[#398]](https://github.com/nf-core/smrnaseq/pull/398) - Update [Input channels](https://github.com/nf-core/smrnaseq/issues/390) - Updated channel and params handling through workflows.
- [[#405]](https://github.com/nf-core/smrnaseq/pull/405) - Fix [Umicollapse algo wrong set](https://github.com/nf-core/smrnaseq/issues/404) - Fix potential bug in Umicollapse (not effective as we do not allow PE data in smrnaseq - but for consistency)
- [[#420]](https://github.com/nf-core/smrnaseq/pull/420) - Fix [mirTrace produces an error in test nextflex](https://github.com/nf-core/smrnaseq/issues/419) - Allow config mode to be used in mirtrace/qc
- [[#425]](https://github.com/nf-core/smrnaseq/pull/425) - Raise [minimum required NXF version for pipeline](https://github.com/nf-core/smrnaseq/issues/424) - usage of `arity` in some modules now requires this
- [[#426]](https://github.com/nf-core/smrnaseq/pull/426) - Add [nf-core mirtop](https://github.com/nf-core/smrnaseq/issues/426) - replace local for nf-core `mirtop`
- [[#427]](https://github.com/nf-core/smrnaseq/pull/427) - Add [nf-core pigz uncompress](https://github.com/nf-core/smrnaseq/issues/422) - replace local `mirdeep_pigz`
- [[#429]](https://github.com/nf-core/smrnaseq/pull/429) - Make [saving of intermediate files optional](https://github.com/nf-core/smrnaseq/issues/424) - Allows user to choose whether to save intermediate files or not. Replaces several params that referred to the same such as `params.save_aligned` and `params.save_aligned_mirna_quant`.
- [[#430]](https://github.com/nf-core/smrnaseq/pull/430) - Emit a [warning if paired-end end data is used](https://github.com/nf-core/smrnaseq/issues/423) - pipeline handles SE data
- [[#432]](https://github.com/nf-core/smrnaseq/pull/432) - Update [MultiQC and all modules to latest version](https://github.com/nf-core/smrnaseq/issues/428) - Include UMIcollapse module in MultiQC.
- [[#435]](https://github.com/nf-core/smrnaseq/pull/435) - Replace local instances of bowtie for nf-core [`bowtie2`](https://github.com/nf-core/smrnaseq/issues/434) and [`bowtie1`](https://github.com/nf-core/smrnaseq/issues/433) - Additionally adds a `bioawk` module that cleans fasta files.
- [[#438]](https://github.com/nf-core/smrnaseq/pull/438) - Update [Mirtop to latest version](https://github.com/nf-core/smrnaseq/issues/437) - Process samples separately and join results with `CSVTK_JOIN`.
- [[#439]](https://github.com/nf-core/smrnaseq/pull/439) - Fix [Fix paired end samples processing](https://github.com/nf-core/smrnaseq/issues/415) - Fix paired end sample handling and add test profile.
- [[#441]](https://github.com/nf-core/smrnaseq/pull/441) - Migrate [local contaminant bowtie to nf-core](https://github.com/nf-core/smrnaseq/issues/436) - Replace local processes with `BOWTIE2_ALIGN`.
- [[#443]](https://github.com/nf-core/smrnaseq/pull/443) - Migrate [mirna and genome_quant bowtie to nf-core](https://github.com/nf-core/smrnaseq/issues/436) - Replace local processes with `BOWTIE_ALIGN`.

## v2.3.1 - 2024-04-18 - Gray Zinc Dalmation Patch

- [[#328]](https://github.com/nf-core/smrnaseq/pull/328) - Fix [casting issue](https://github.com/nf-core/smrnaseq/issues/327) in mirtrace module
- [[#334]](https://github.com/nf-core/smrnaseq/pull/334) - Fix [input channel cardinality](https://github.com/nf-core/smrnaseq/issues/331) in `MIRDEEP2_RUN` module
- [[#334]](https://github.com/nf-core/smrnaseq/pull/334) - Fix [bowtie conda version](https://github.com/nf-core/smrnaseq/issues/333) in `BOWTIE_MAP_SEQ` module
- [[#335]](https://github.com/nf-core/smrnaseq/pull/335) - Final fix for [casting issue](https://github.com/nf-core/smrnaseq/issues/327) in mirtrace module
- [[#337]](https://github.com/nf-core/smrnaseq/pull/337) - Fix [three_prime_adapter issue](https://github.com/nf-core/smrnaseq/issues/326), allow `auto-detect` as value
- [[#342]](https://github.com/nf-core/smrnaseq/pull/342) - Fix [phred offset issue](https://github.com/nf-core/smrnaseq/issues/341), allow specifying phred offset for FASTQ files
- [[#343]](https://github.com/nf-core/smrnaseq/pull/343) - Fix [mirdeep2 output missing](https://github.com/nf-core/smrnaseq/issues/330), fix mirdeep2 outputs missing in outdir

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `multiqc`  | 1.20        | 1.21        |

## v2.3.0 - 2024-02-23 - Gray Zinc Dalmatian

- [[#307]](https://github.com/nf-core/smrnaseq/pull/307) - Clean up config file and improve output folder structure
- [[#299]](https://github.com/nf-core/smrnaseq/issues/299) - Bugfix for missing inputs in BAM stats (`genome_quant.r`)
- [[#164]](https://github.com/nf-core/smrnaseq/pull/164) - UMI Handling Feature implemented in the pipeline
- [[#302]](https://github.com/nf-core/smrnaseq/pull/302) - Merged in nf-core template v2.11.1
- [[#294]](https://github.com/nf-core/smrnaseq/pull/294) - Fixed contamination screening issues
- [[#309]](https://github.com/nf-core/smrnaseq/pull/309) - Merged in nf-core template v2.12.0
- [[#310]](https://github.com/nf-core/smrnaseq/pull/310) - Removed unnecessarily separate mirtrace subworkflow, now using module instead
- [[#311]](https://github.com/nf-core/smrnaseq/pull/311) - Fix use of FASTP, set `three_prime_adapter` per default
- [[#314]](https://github.com/nf-core/smrnaseq/pull/314) - Add parameters to control publishing of intermediate results
- [[#317]](https://github.com/nf-core/smrnaseq/pull/317) - Fixed issue with bowtie indices directly supplied
- [[#318]](https://github.com/nf-core/smrnaseq/pull/318) - Merged in nf-core template v2.13.0 and pinned nf-validator

### Parameters

| Old parameter | New parameter                |
| ------------- | ---------------------------- |
|               | `--with_umi`                 |
|               | `--umitools_extract_method`  |
|               | `--umitools_method`          |
|               | `--skip_umi_extract`         |
|               | `--umitools_bc_pattern`      |
|               | `--umi_discard_read`         |
|               | `--save_umi_intermeds`       |
|               | `--save_aligned`             |
|               | `--save_aligned_mirna_quant` |
|               | `--save_merged`              |

### Software dependencies

| Dependency    | Old version | New version |
| ------------- | ----------- | ----------- |
| `multiqc`     | 1.15        | 1.20        |
| `edgeR`       | 3.36.0      | 4.0.2       |
| `limma`       | 3.50.0      | 3.58.1      |
| `bioconvert`  | 0.4.3       | 1.1.1       |
| `mirdeep`     | 2.0.1       | 2.0.1.3     |
| `seqkit`      | 2.3.1       | 2.6.1       |
| `fastqc`      | 0.11.4      | 0.12.1      |
| `samtools`    | 1.17        | 1.19.2      |
| `umitools`    | <none>      | 1.1.5       |
| `umicollapse` | <none>      | 1.0.0       |

## [v2.2.4](https://github.com/nf-core/smrnaseq/releases/tag/2.2.4) - 2023-11-03

- Update template to 2.10
- [[#289]](https://github.com/nf-core/smrnaseq/issues/289) - Bugfix for issue with mirdeep2 channels ()
- [[#288]](https://github.com/nf-core/smrnaseq/issues/288) - Bugfix for issue with handling malformed GFF3 from mirbase
- Updated dependencies, including FASTQC, MultiQC 1.17, fastP and samtools to latest versions

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

## [v2.2.3](https://github.com/nf-core/smrnaseq/releases/tag/2.2.3) - 2023-09-06

- [[#271]](https://github.com/nf-core/smrnaseq/issues/271) - Bugfix for parsing hairpin and mature fasta files

## [v2.2.2](https://github.com/nf-core/smrnaseq/releases/tag/2.2.2) - 2023-09-04

- [[#253]](https://github.com/nf-core/smrnaseq/pull/253) - Remove globs from process alias when using ECR containers
- [[#237]](https://github.com/nf-core/smrnaseq/issues/237) - Fix illumina protocol clip parameters to default
- Remove public_aws_ecr profile
- [[#269]](https://github.com/nf-core/smrnaseq/pull/269) - Updated miRBase URLs with new location (old ones were broken)

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `multiqc`  | 1.13        | 1.15        |
| `fastp`    | 0.23.2      | 0.23.4      |

## [v2.2.1](https://github.com/nf-core/smrnaseq/releases/tag/2.2.1) - 2023-05-12

### Enhancements & fixes

- [[#238](https://github.com/nf-core/smrnaseq/issues/238)] - Restored the missing mirtop outputs
- [[#242](https://github.com/nf-core/smrnaseq/issues/242)] - Fixed mirtrace using wrong input fastq files (untrimmed)

## [v2.2.0](https://github.com/nf-core/smrnaseq/releases/tag/2.2.0) - 2023-04-26

- [[#220](https://github.com/nf-core/smrnaseq/issues/220)] - Fixed an issue where miRTrace reports fastq basename instead of sample ID
- [[#208](https://github.com/nf-core/smrnaseq/issues/208)] - Fixed usability of `--skip_qc` parameter
- Updated FASTP module to allow direct addition of adapter_fasta file to it
- [[#205](https://github.com/nf-core/smrnaseq/issues/205)] - Fix mirTrace Report to be a single report again instead of per sample reports
- [[#206](https://github.com/nf-core/smrnaseq/issues/206)] - Use % as separator in sed commands to enable conda working properly on OSX and Linux
- [[#207](https://github.com/nf-core/smrnaseq/issues/224)] - Fix Samplesheet print error
- Group samples by adapter sequence before running mirtrace
- Remove `--skip_qc` parameter

## [v2.1.0](https://github.com/nf-core/smrnaseq/releases/tag/2.1.0) - 2022-10-20 Maroon Tin Dalmatian

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
- [#49](https://github.com/nf-core/smrnaseq/issues/49) - Integrated the existing umitools modules into the pipeline and extend the deduplication step.
- [#55](https://github.com/lpantano/seqcluster/pull/55) - update seqcluster to fix UMI-detecting bug

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

### Parameters

| Old parameter        | New parameter    |
| -------------------- | ---------------- |
| `--conda`            | `--enable_conda` |
| `--clusterOptions`   |                  |
| `--publish_dir_mode` |                  |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

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
