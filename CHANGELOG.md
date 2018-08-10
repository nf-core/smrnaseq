# nf-core/smrnaseq

## [1.0](https://github.com/nf-core/smrnaseq/releases/tag/1.0) - 2018-08-06
* Switch from SciLifeLab/NGI-smRNAseq to nf-core/smrnaseq.
* Use Bowtie 1 instead of Bowtie 2 for the alignment to host reference genome.
* Add option for sequencing centre in BAM file.
* Software updates: trim-galore 0.4.5 to 0.5.0; samtools 1.8 to 1.9; multiqc 1.5 to 1.6.

Change bowtie parameters from `bowtie -n x -l 15 -k 10 --best (x=0 for mature and x=1 for hairpin)` into `bowtie -k 1 -m 1 --best --strata`.

## [0.1dev](https://github.com/SciLifeLab/NGI-smRNAseq/releases/tag/0.1dev) - 2018-05-14
* Change bowtie parameters from `bowtie -n x -l 15 -k 10 --best (x=0 for mature and x=1 for hairpin)` into `bowtie -k 1 -m 1 --best --strata`.
