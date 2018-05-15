# NGI-smRNAseq

## 2018-05-14
Change bowtie parameters from `bowtie -n x -l 15 -k 10 --best (x=0 for mature and x=1 for hairpin)` into `bowtie -k 1 -m 1 --best --strata`.
