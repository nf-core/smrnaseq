# nf-core/smrnaseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/smrnaseq/usage](https://nf-co.re/smrnaseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

### `protocol`

This option indicates the experimental protocol used for the sample preparation. Currently supporting:

- 'illumina': three_prime_adapter (`TGGAATTCTCGGGTGCCAAGG`), clip_r1 (`0`), three_prime_clip_r1 (`0`)
- 'nextflex': NEXTFLEX Small RNA-Seq Kit v3; three_prime_adapter (`TGGAATTCTCGGGTGCCAAGG`), clip_r1 (`4`), three_prime_clip_r1 (`4`)
- 'nextflex_v4': NEXTFLEX Small RNA-Seq Kit v4; three_prime_adapter (`TGGAATTCTCGGGTGCCAAGG`), clip_r1 (`0`), three_prime_clip_r1 (`0`), fastp_min_length (`16`)
- 'qiaseq': three_prime_adapter (`AACTGTAGGCACCATCAAT`), clip_r1 (`0`), three_prime_clip_r1 (`0`), mirna_forward_only (`true`); with `--with_umi`, also the QIAseq UMI pattern (see [UMI handling](#umi-handling))
- 'cats': three_prime_adapter (`AAAAAAAA`), clip_r1(`3`), three_prime_clip_r1 (`0`)

This option is not chosen as a parameter but as an additional profile that sets the corresponding `three_prime_adapter`, `clip_r1` and `three_prime_clip_r1` parameters and, where applicable, `fastp_min_length` and `mirna_forward_only` accordingly. You can choose to either use any of the provided profiles by running the pipeline with e.g. `illumina` to set the defaults as described above in a more convenient way.

```bash
-profile your_other_profiles,illumina
```

For NEXTFLEX libraries, use `nextflex` for the NEXTFLEX Small RNA-Seq Kit v3 chemistry and `nextflex_v4` for the v4 chemistry. For v4, the settings follow the [manufacturer's trimming instructions](https://resources.revvity.com/pdfs/man-nextflex-small-rna-trimming-instructions.pdf), which specify the 3' adapter and a minimum read length of 16 nt without additional fixed 5' or 3' clipping.

In case you have a custom protocol, please supply the `three_prime_adapter`, `clip_r1` and `three_prime_clip_r1` manually and adapt `fastp_min_length` if required by your protocol.

The parameter `--three_prime_adapter` is set to the Illumina TruSeq single index adapter sequence `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`. This is also to ensure, that the auto-detect functionality of `FASTP` is disabled. Please make sure to adapt this adapter sequence accordingly for your run.

:warning: If you do not choose a profile that sets the `three_prime_adapter`, `clip_r1` and `three_prime_clip_r1` options, the pipeline won't run. If you want to auto-detect the adapters using `fastp`, please set `--three_prime_adapter` to `auto-detect`.

miRTrace cannot auto-detect adapters. With `--three_prime_adapter auto-detect`, miRTrace QC runs on the raw reads without adapter trimming, so its read length, RNA type and contamination statistics are computed on untrimmed reads. To get adapter-trimmed miRTrace QC, set `--three_prime_adapter` to the adapter sequence of your library. As fastp's adapter auto-detection is also unreliable on small RNA data, setting `--three_prime_adapter` explicitly to your kit's adapter is recommended for meaningful miRTrace QC.

### `mirtrace_species` or `mirgenedb_species`

It should point to the 3-letter species name used by [miRBase](https://www.mirbase.org/browse) or [MirGeneDB](https://www.mirgenedb.org/browse). Note the difference in case for the two databases.

`mirtrace_species` is required when using miRBase (the default). Besides selecting the miRTrace species, it is used to select the miRNAs of your species from the `mature` and `hairpin` FASTA files and, unless `mirna_gtf` is given, to download the miRBase GFF3 file for mirtop (see below). Use the miRBase species code, e.g. `hsa`, `mmu` or `oar` for sheep. It cannot be set to a placeholder such as `none` to skip miRTrace, because the miRNA quantification would then find no miRNAs.

### miRNA related files

Different parameters can be set for the two supported databases. By default `miRBase` will be used with the parameters below.

- `mirna_gtf`: If not supplied by the user, then `mirna_gtf` will point to the latest GFF3 file in miRBase: `https://mirbase.org/download/${params.mirtrace_species}.gff3`
- `mature`: points to the FASTA file of mature miRNA sequences. Default: `https://mirbase.org/download/mature.fa`
- `hairpin`: points to the FASTA file of precursor miRNA sequences. Default: `https://mirbase.org/download/hairpin.fa`

#### Species without a miRBase GFF3 file

miRBase only provides GFF3 files for species whose miRNAs are mapped to a genome assembly. For many species in miRBase, such as sheep (`oar`) or rabbit (`ocu`), the file `https://mirbase.org/download/${params.mirtrace_species}.gff3` does not exist and the pipeline stops with an error that the miRBase GFF3 file could not be found. These species can still be analysed: their miRNAs are in the miRBase `mature` and `hairpin` FASTA files, and miRTrace supports them (`mirtrace --list-species` lists the supported species codes).

The GFF3 file is only used by mirtop for the isomiR annotation (`mirna_quant/mirtop/`). You can either:

- skip mirtop with `--skip_mirtop`. No GFF3 file is required or downloaded, and miRTrace QC, the mature and hairpin miRNA quantification (`mirna_quant/edger_qc/`), the genome alignment and miRDeep2 still run. For example, for sheep:

  ```bash
  nextflow run nf-core/smrnaseq \
      -profile docker,illumina \
      --input samplesheet.csv \
      --outdir results \
      --mirtrace_species oar \
      --fasta Ovis_aries.ARS-UI_Ramb_v2.0.dna.toplevel.fa.gz \
      --skip_mirtop
  ```

- or provide a GFF3 file with `--mirna_gtf`. It must be in the miRBase GFF3 format (`miRNA_primary_transcript` and `miRNA` features with `ID`, `Name` and `Derives_from` attributes, and `Name` values matching the miRNA names in the `hairpin` and `mature` FASTA files), as written by miRBase for other species. A gene annotation GTF, e.g. from Ensembl, does not work.

If your species is in [MirGeneDB](https://mirgenedb.org/browse), you can also use MirGeneDB, which provides a GFF file for every species it contains (see below).

Without `--genome`, set `--fasta` (or `--bowtie_index`) to a genome of your species for the genome alignment and miRDeep2. The `mature` and `hairpin` defaults already contain all miRBase species, so they do not need to be changed.

If MirGeneDB should be used instead it needs to be specified using `--mirgenedb` and use the parameters below.

- `mirgenedb_gff`: The GFF file cannot be downloaded automatically due to the presence of short-term tokens in the URLs. Therefore, the user must manually provide the GFF file, either for their species of interest or for all species, by downloading it from [MirGeneDB](https://mirgenedb.org/download). The provided dataset will be automatically filtered based on the species specified with the `--mirgenedb_species` parameter.
- `mirgenedb_mature`: This parameter should point to the FASTA file containing mature miRNA sequences. The file can be manually downloaded from [MirGeneDB](https://mirgenedb.org/download).
- `mirgenedb_hairpin`: This parameter should point to the FASTA file containing precursor miRNA sequences. Note that MirGeneDB does not offer a dedicated hairpin file, but the precursor sequences can be downloaded from [MirGeneDB](https://mirgenedb.org/download) and used instead.

See examples in: [the test-datasets repository of nf-core](https://github.com/nf-core/test-datasets/tree/smrnaseq/MirGeneDB).

> [!NOTE]
> `mirtop` reads the name of the miRNA database from the header of the GFF3 file, so the GFF3 passed via `--mirna_gtf` must be a plain-text file that starts with the miRBase header comments (e.g. `# microRNAs: miRBase v22.1`). Otherwise `MIRTOP_GFF` fails with `Database not found in --mirna <species>.gff3. Use --database argument to add a custom source.`
>
> - When downloading the file yourself, use `https://mirbase.org/download/<species>.gff3` (e.g. `wget https://mirbase.org/download/mmu.gff3`). The `https://mirbase.org/download/CURRENT/genomes/<species>.gff3` URL serves the file wrapped in HTML (the first line looks like `<p>##gff-version 3<br>...`), which `mirtop` cannot parse.
> - For a custom miRBase-style GFF3 without these header lines, either add a comment line such as `# microRNAs: miRBase v22.1` to the top of the file, or pass the database name to `mirtop` with a custom config (`-c custom.config`):
>
> ```groovy
> process {
>     withName: '.*:BAM_STATS_MIRNA_MIRTOP:MIRTOP_GFF' {
>         ext.args = '--database miRBase'
>     }
> }
> ```
>
> The database name must contain `miRBase` for coordinates of `miRNA_primary_transcript` features to be recognised. With `--mirgenedb`, the pipeline already passes `--database MirGeneDB` automatically.

> [!NOTE] > `mirtop` is hard-coded to use the `pre` sequences, which originate from the hairpin FASTA, rather than the `pri` sequences, which come from the mature FASTA. Users must provide `pre` files from the start to ensure consistency between the FASTA and GFF files, as the coordinates in the GFF file are referenced to `pre` sequences. This also ensures that names in the BAM file will match those in the GFF.

### miRNA reference orientation

The mature and hairpin FASTA files contain miRNA sequences in RNA (sense) orientation. In all library protocols that the pipeline has a profile for, read 1 starts at the 5′ end of the small RNA, so after trimming a read has the same orientation as the miRNA it comes from. By default, the pipeline still lets Bowtie 1 align reads to both strands of the mature, hairpin and seqcluster-to-hairpin references. A read from an abundant miRNA can then also align to the reverse complement of another annotated miRNA. Several miRBase entries are antisense to abundant miRNAs, for example miR-3184-3p (miR-423-5p), miR-103b (miR-103a-3p and miR-107), miR-3529-3p (miR-7-5p), miR-122b-3p (miR-122-5p) and miR-451b (miR-451a). Some precursors (e.g. mir-486) are self-complementary, so a read from one arm also aligns to the other arm in reverse.

`--mirna_forward_only` adds Bowtie 1 [`--norc`](https://bowtie-bio.sourceforge.net/manual.shtml#the-nofw-norc-options) to these three alignments, so only forward alignments to the miRNA references are reported. The mismatch and multimapping settings are unchanged. Reads that no longer align to the mature miRNAs continue to the hairpin and genome alignments. Genome, miRDeep2 and contaminant alignments are not changed: this is about the orientation relative to the RNA reference, not the genomic strand, and a miRNA encoded on the minus strand of the genome still aligns forward to its own sequence.

The effect on the outputs:

- `mirna_quant/edger_qc/` (`mature_counts.csv`, `hairpin_counts.csv` and the files derived from them): these counts come from `samtools idxstats`, which counts every alignment record, including reverse-strand ones. With `--mirna_forward_only`, the antisense miRNAs listed above are no longer reported with counts that actually come from their sense partners. In the pipeline test data (`illumina`, `nextflex` and `qiaseq` with UMIs), mature totals dropped by 1.7-2.9 % and hairpin totals by 2-31 %, and 16-74 miRNAs or hairpins per sample that only had reverse-strand alignments were no longer reported. The most abundant miRNAs did not change.
- `mirna_quant/mirtop/` (`mirna.tsv`, `joined_samples_mirtop.tsv`): mirtop already ignores reverse-strand alignments, so these counts change very little (below 0.1 % in the test data). A few miRNAs can gain counts, because a reverse-strand alignment with fewer mismatches no longer hides a forward alignment of the same read (`--best --strata`).

The `qiaseq` profile enables `--mirna_forward_only`, in line with QIAGEN's [strand-specific QIAseq miRNA quantification](https://resources.qiagenbioinformatics.com/manuals/rnaanalysisportal/4100/index.php?manual=QIAseq_miRNA_Library_Kit.html) (QIAGEN's [miRNA quantification](https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/2104/index.php?manual=Quantify_miRNA.html) only searches the plus strand of the RNA reference when strand-specific). This also applies with `--with_umi`: UMI extraction moves the UMI to the read name and keeps the insert in its original orientation.

For the other protocols, the option is not enabled by default (to keep results comparable with earlier releases), but it is safe and recommended to set `--mirna_forward_only` for:

| Profile                   | Library chemistry                                              | Read 1 orientation               |
| ------------------------- | -------------------------------------------------------------- | -------------------------------- |
| `illumina`                | TruSeq Small RNA: 5′ and 3′ adapter ligation                   | sense                            |
| `nextflex`, `nextflex_v4` | NEXTFLEX Small RNA-Seq v3/v4: 5′ and 3′ adapter ligation       | sense                            |
| `qiaseq`                  | QIAseq miRNA: 3′ and 5′ adapter ligation, UMI in the RT primer | sense (enabled by default)       |
| `cats`                    | CATS: 3′ polyadenylation and template switching at the 5′ end  | sense (`GGG` + insert + poly(A)) |

Do not use `--mirna_forward_only` for reads in antisense orientation, for example if read 2 of a paired-end run was provided as `fastq_1`, or for a custom protocol with unknown or mixed orientation. In these cases, a large fraction of reads aligning in reverse to the mature miRNAs without the option (see `mirna_quant/bam/mature/` with `--save_intermediates`) indicates the problem. For paired-end samplesheets, the pipeline only uses `fastq_1`.

> [!WARNING]
> With the `qiaseq` profile, the edgeR count tables differ from earlier releases (see above). To reproduce the previous both-strand miRNA alignment, set `--mirna_forward_only false`. Rerun the alignment rather than filtering reverse-strand records from existing BAM files, because restricting the search can reveal forward alignments that were hidden before.

### Genome

- `fasta`: the reference genome FASTA file
- `bowtie_index`: points to the folder containing the `bowtie` indices for the genome reference specified by `fasta`.
- `bowtie_align_args`: Optional flags to pass to Bowtie 1 for the genome mapping, e.g. `--bowtie_align_args '-v 1 --best'`. Only a subset of Bowtie 1 flags are allowed to keep compatibility with the pipeline. If invalid or disallowed values are included, the pipeline will exit with an error message. A value holding a single flag needs the `=` form, for example `--bowtie_align_args=--best`, because `--bowtie_align_args --best` is parsed as the boolean `true` and is rejected. Short flags take their value as a separate word (`-v 1`, not `-v=1`).
  - These flags are appended after the pipeline's defaults for genome mapping (`-t -k 50 --best --strata -e 99999 --chunkmbs 2048`). When a flag is repeated, Bowtie uses the last value, so for example `--bowtie_align_args '-k 10'` replaces `-k 50`.
  - `-k` must be at least 2, because the defaults include `--strata`, which Bowtie only accepts with `-k` above 1. To report a single alignment per read, override `ext.args` for `NFCORE_SMRNASEQ:GENOME_QUANT:BOWTIE_MAP_GENOME` instead (see [Custom Tool Arguments](#custom-tool-arguments)).

> [!NOTE]
> if the FASTA file in `fasta` is not the same file used to generate the `bowtie` indices, then the pipeline will fail.

If no `bowtie_index` is given, the pipeline builds one from `fasta`. Add `--save_reference` to save it, together with the cleaned FASTA it was built from, to `<outdir>/genome/` and reuse the index in future runs with `--bowtie_index <outdir>/genome/index/bowtie`, together with the same `--fasta` (or the saved `_clean.fa`). `--save_reference` also saves the other indices the pipeline builds (mature and hairpin miRNA Bowtie indices and, with `--filter_contamination`, the contaminant Bowtie2 indices); see [the output documentation](https://nf-co.re/smrnaseq/output#reference-genome-files).

#### Genome quantification and miRDeep2

Genome quantification and miRDeep2 are optional steps. There is no separate parameter to switch genome quantification on: it runs whenever the pipeline has a genome FASTA, which you can provide in either of two ways:

- `--fasta /path/to/genome.fa`, or
- `--genome <ID>` with an iGenomes reference (for example `--genome GRCh38`), which sets `fasta` from [`conf/igenomes.config`](../conf/igenomes.config). An explicit `--fasta` takes precedence.

Without a genome FASTA, both steps are skipped without an error, and the pipeline prints a warning at start-up (which only mentions genome quantification when `--skip_mirdeep` is set). This also applies when only `--bowtie_index` is given: the genome Bowtie index is used only together with `fasta`, and is otherwise ignored, which the warning then points out. With `--igenomes_ignore`, `--genome` does not set a FASTA, so you have to pass `--fasta`.

The genome Bowtie 1 index comes from `--bowtie_index` (a directory or a `.tar.gz` archive) if given, and is otherwise built from the FASTA. The iGenomes references used by this pipeline do not provide a Bowtie 1 index, so with `--genome` the index is built at the start of every run, which takes time for large genomes. Build it once with `--save_reference` and pass it back with `--bowtie_index` in later runs, as described above.

When a genome FASTA is available, the pipeline runs:

1. **Genome quantification** (`GENOME_QUANT`): the reads that did not align to the mature miRNA or hairpin references are aligned to the genome with Bowtie 1 (see `bowtie_align_args` above for the default and extra options), then sorted, indexed and summarised with SAMtools (`stats`, `flagstat` and `idxstats`). The alignment statistics are added to the MultiQC report. This step does not produce a gene- or feature-level count table: `idxstats` reports the number of reads per genome sequence (chromosome or contig). miRNA counts come from the [miRNA quantification](https://nf-co.re/smrnaseq/output#edger) steps.
2. **miRDeep2** for known and novel miRNA discovery, using the same genome FASTA and Bowtie 1 index. miRDeep2 runs on all trimmed (and, with `--filter_contamination`, filtered) reads, not only on those left after miRNA alignment. Add `--skip_mirdeep` to run genome quantification without miRDeep2.

For example:

```bash
nextflow run nf-core/smrnaseq \
    -profile docker,illumina \
    --input samplesheet.csv \
    --outdir results \
    --mirtrace_species hsa \
    --fasta /path/to/genome.fa \
    --bowtie_index /path/to/bowtie_index/
```

The SAMtools statistics are always saved to `<outdir>/genome_quant/bam/`. The genome alignments (BAM files), the Bowtie logs and the reads that did not align to the genome are saved there only with `--save_intermediates`. See [the output documentation](https://nf-co.re/smrnaseq/output#bowtie) for the full list of files, and [miRDeep2](https://nf-co.re/smrnaseq/output#mirdeep2) for its results.

### Contamination filtering

This step has, until now, only been tested for human data. Unexpected behaviour can occur when using it with a different species.

Contamination filtering of the sequencing reads is optional and can be invoked using the `filter_contamination` parameter. FASTA files with

- `rrna`: Used to supply a FASTA file containing rRNA contamination sequence.
- `trna`: Used to supply a FASTA file containing tRNA contamination sequence. e.g. `http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa`
- `cdna`: Used to supply a FASTA file containing cDNA contamination sequence. e.g. `ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz` The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `ncrna`: Used to supply a FASTA file containing ncRNA contamination sequence. e.g. `ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz` The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `pirna`: Used to supply a FASTA file containing piRNA contamination sequence. e.g. The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `other_contamination`: Used to supply an additional filtering set. The FASTA file is first compared to the available miRNA sequences and overlaps are removed.

## mirDeep2

If the software encounters an error with exit status 255, it will be ignored, and the pipeline will continue to complete. In such cases, the pipeline will log a note that includes the path to the work directory where the issue occurred. You can inspect this work directory to examine your input data and troubleshoot the issue.

Error 255 is typically related to the core algorithm of miRDeep generating empty output files. This often happens when the reads being processed do not correspond to putative mature miRNA sequences, or if the provided precursors do not meet the criteria for valid miRNA precursors, both of which may stem from the input reads used. A common cause of this error is running the pipeline with a small subset of the input reads.

### UMI handling

With `--with_umi`, the pipeline removes PCR duplicates using unique molecular identifiers (UMIs) before any alignment. The steps run in this order:

1. `umi_tools extract` runs on the raw reads. It uses `--umitools_extract_method` and `--umitools_bc_pattern` to move the UMI into the read name and to remove the sequences around it, so that only the insert of the RNA molecule remains. Reads that do not match the pattern are discarded and counted in the `umi_tools extract` log.
2. `fastp` filters the inserts by quality and length (`--fastp_min_length`, `--fastp_max_length`) and clips `--clip_r1` bases from the 5' end; `--three_prime_clip_r1` bases are then clipped from the 3' end. With `--umitools_extract_method regex`, the pattern is expected to remove the kit adapter, so `fastp` does not trim adapters. With `string`, `fastp` trims `--three_prime_adapter` as without UMIs. The pipeline's list of known miRNA adapters is not used with UMIs.
3. The reads are deduplicated on **UMI + insert sequence** with [UMICollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse) in FastQ mode. UMICollapse cannot read the UMI from the read name, so the pipeline moves the UMI to the start of the read sequence, runs UMICollapse, and removes the UMI from the sequence again. By default, reads whose UMI + insert sequences differ by one substitution are merged (see [Choosing a UMI strategy](#choosing-a-umi-strategy)). Reads of different lengths are never merged.

The deduplicated inserts are then used for all later steps; the UMI stays in the read name. miRTrace QC runs on the raw reads as without UMIs.

> [!NOTE]
> Before nf-core/smrnaseq 2.5.0, UMICollapse ran on the adapter-trimmed reads before UMI extraction, so it compared insert, kit adapter and UMI together, and `umi_tools extract` ran afterwards. Extracting the UMI first removes sequencing errors in the kit adapter and in sequence after the UMI from the comparison ([#560](https://github.com/nf-core/smrnaseq/issues/560)). `-profile qiaseq` now also works with `--with_umi`.

#### QIAseq miRNA

QIAseq miRNA reads have the structure `<insert><AACTGTAGGCACCATCAAT><12 nt UMI><Illumina adapter>`. `-profile qiaseq` sets the matching UMI pattern, so `--with_umi` is all that is needed:

```bash
nextflow run nf-core/smrnaseq \
    -profile docker,qiaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --mirtrace_species hsa \
    --with_umi
```

The profile sets `--umitools_extract_method regex` and `--umitools_bc_pattern '.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)'`: the UMI is the 12 nt after the kit adapter (allowing up to two mismatches in the adapter), and the adapter, the UMI and everything after it are removed. The profile's `--three_prime_adapter` (the kit adapter) is used by miRTrace on the raw reads.

#### Other UMI read structures

For other kits, set `--umitools_extract_method` and `--umitools_bc_pattern` so that only the insert remains after extraction (see the [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#extract-method) and your kit's documentation for its read structure):

- **UMI at a fixed position at the 5' end** of the read: for example, `--umitools_extract_method string --umitools_bc_pattern NNNNNNNN` extracts an 8 nt UMI. Set `--three_prime_adapter` to your kit's 3' adapter: `fastp` removes it after extraction.
- **UMI next to an adapter or other kit sequence** (as for QIAseq): use `--umitools_extract_method regex` with a pattern that finds the UMI relative to the kit sequence and discards everything after the insert, as in the QIAseq pattern above. With `regex`, `fastp` does not trim adapters, so anything the pattern leaves in the read stays there.

`--clip_r1` and `--three_prime_clip_r1` are applied after UMI extraction, so they count bases of the insert, not of the UMI.

#### UMI parameters

- `--with_umi`: enable UMI deduplication. Requires `--umitools_bc_pattern` (set by `-profile qiaseq`).
- `--umitools_extract_method`: `string` (default) or `regex`, passed to `umi_tools extract --extract-method`.
- `--umitools_bc_pattern`: the UMI pattern, passed to `umi_tools extract --bc-pattern`.
- `--umitools_method`: despite its name, this sets the UMICollapse deduplication algorithm (`--algo`): `dir` (directional, default), `adj` (adjacency) or `cc` (connected components). All three merge UMIs that differ by one substitution. See [UMI deduplication algorithm](#umi-deduplication-algorithm).
- `--save_umi_intermeds`: save the reads after `umi_tools extract` and after deduplication (see the [output documentation](https://nf-co.re/smrnaseq/output#umi-tools-extract)).
- `--skip_umi_extract_before_dedup` and `--umi_discard_read`: deprecated and ignored. The UMI is always extracted before deduplication, and only read 1 is used. The pipeline prints a warning if they are set.

#### UMI deduplication algorithm

UMICollapse groups reads whose UMI + insert sequences differ by at most one substitution (by default) and keeps one read per group. How the groups are formed is set by `--umitools_method`. Despite its name, this parameter is passed to UMICollapse as `--algo`; it does not affect UMI-tools. The available algorithms match the UMI-tools grouping methods of the same name:

- `dir` (default): directional. Two sequences are grouped when the more abundant one has at least twice the count of the other (minus one). This is the most conservative choice against merging distinct molecules.
- `adj`: adjacency. Each cluster of connected sequences is resolved into groups around the most abundant sequences.
- `cc`: connected components. All sequences connected by single mismatches form one group. This collapses the most reads.

For example:

```bash
--with_umi --umitools_method adj
```

Other UMICollapse options, such as the mismatch threshold (`-k`, default 1) or the directional threshold (`-p`, default 0.5), are not exposed as pipeline parameters. You can set them in a custom config file passed with `-c`. Note that `ext.args` replaces the arguments set by the pipeline, so keep the `--algo` option. For example, to allow two mismatches (for exact UMI matching with `-k 0`, see [Choosing a UMI strategy](#choosing-a-umi-strategy)):

```groovy
process {
    withName: '.*:UMICOLLAPSE_FASTQ' {
        ext.args = '--algo dir -k 2'
    }
}
```

See the [UMICollapse documentation](https://github.com/Daniel-Liu-c0deb0t/UMICollapse#command-line-arguments) for all options.

### Choosing a UMI strategy

UMI deduplication has to decide when two reads come from the same molecule. Reads with the same insert and UMI are duplicates. Reads whose UMIs differ by one base can be duplicates with a PCR or sequencing error in the UMI, or two different molecules whose UMIs happen to be similar. The pipeline offers two strategies:

- **One-mismatch clustering (default):** UMICollapse `--algo dir` (directional) merges a read into a more abundant read whose UMI + insert differs by one substitution. This corrects UMI errors, but once a single insert sequence has been seen with a large part of all possible UMIs, distinct molecules have neighbouring UMIs by chance and are merged as well (_saturation_). The counts of the most abundant miRNAs are then underestimated, and more so the deeper the library is sequenced.
- **Exact UMI matching:** only reads with identical UMI + insert are merged. This is not affected by saturation, but reads with a UMI error are counted as separate molecules: a near-uniform overcount of a few percent, which mostly cancels in library-size normalisation. Set it with a custom config (see [Custom Tool Arguments](#custom-tool-arguments)):

  ```groovy
  process {
      withName: '.*:UMICOLLAPSE_FASTQ' {
          ext.args = '--algo dir -k 0'
      }
  }
  ```

  Keep the `--algo` option: `ext.args` replaces the pipeline's arguments. Exact matching is also much faster (see below).

**Benchmark.** The strategies were compared on simulated QIAseq-like libraries with known molecules (500,000–1,000,000 molecules with 40% duplicates, 1% UMI and insert errors per PCR copy, 0.2% sequencing errors per base, and one miRNA with 40–60% of all molecules), and on the `test_umi` sample (1.1M reads) and a 4.5M-read subsample of the same QIAseq plasma library:

|                                                                         | One mismatch (`dir`, default) | Exact UMI (`-k 0`)                   | Previous order (UMICollapse before extraction, `dir`) |
| ----------------------------------------------------------------------- | ----------------------------- | ------------------------------------ | ----------------------------------------------------- |
| Count error per miRNA, simulated, all but the most abundant miRNA       | 0.6–1.1%                      | +8–9% (overcount)                    | 1.1–1.6%                                              |
| Most abundant simulated miRNA (200,000 / 600,000 molecules)             | −7% / −11%                    | +8% / +8%                            | −6% / −10%                                            |
| Real data: most abundant miRNA (miR-16-5p), exact UMI relative to `dir` |                               | 1.9× (1.1M reads), 2.2× (4.5M reads) |                                                       |
| Real data: other miRNAs, relative to the previous order                 | −2% to +1%                    | +4% to +13%                          |                                                       |
| CPU time and peak memory, 4.5M reads, extraction to deduplicated reads  | 720 s, 1.6 GB                 | 306 s, 1.9 GB                        | 761 s, 2.5 GB                                         |
| CPU time of the deduplication step alone, 4.5M reads                    | 436 s                         | 21 s                                 | 669 s                                                 |

How to choose:

- The default (`dir`) is the most accurate for typical miRNAs as long as no insert sequence uses more than about **0.5% of all possible UMIs** (4^UMI length; 84,000 molecules for a 12 nt UMI).
- The pipeline checks this after deduplication: if the most frequent insert sequence of a sample has more deduplicated reads than 0.5% of the UMI space, it prints a **UMI saturation warning**. The statistics for each sample are saved in `umi_dedup/saturation/` (see the [output documentation](https://nf-co.re/smrnaseq/output#umi-collapse-deduplicate)). The warning appears for both `test_umi` samples (1.1M reads each): after deduplication, the canonical miR-16-5p sequence still has about 172,000 reads, 1.0% of the 12 nt UMI space.
- If you see the warning and the abundance of the most abundant miRNAs matters for your analysis (for example for differential expression), consider exact UMI matching. Its overcount is roughly uniform across miRNAs, whereas the undercount from saturation affects only the most abundant ones and grows with depth.

**Runtime and memory.** UMICollapse holds all reads of a sample in memory and its one-mismatch search slows down at saturated sequences, so it is the slowest step of the UMI path; give `UMICOLLAPSE_FASTQ` more memory for deep libraries. `umi_tools extract` is single-threaded (about 270 s CPU per 4.5M reads).

**Planned.** A position-based deduplication backend using [fgumi](https://github.com/fulcrumgenomics/fgumi) after alignment, with [chelae](https://github.com/fulcrumgenomics/chelae) as an optional fast UMI extractor, is planned for nf-core/smrnaseq 2.6 ([#535](https://github.com/nf-core/smrnaseq/issues/535), [#530](https://github.com/nf-core/smrnaseq/issues/530)).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 2 columns ("sample" and "fastq_1"), and a header row as shown in the examples below.

If a second fastq file is provided using another column, the extra data are ignored by this pipeline. The smRNA species should be sufficiently contained in the first read, and so the second read is superfluous data in this smRNA context.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers should match between runs of resequenced samples. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet must have at least 2 columns (`sample` and `fastq1`). A third column can be added if the sample is paired-end (`fastq2`).

> [!NOTE]
> Most of the tools used can't accommodate paired end reads, so whenever paired-end samples are used as inputs, only the R1 files are used by the pipeline.

A final samplesheet file consisting of single-end data and may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz
```

| Column    | Description                                                                                                                                                                            | Requirement |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). | Mandatory   |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             | Mandatory   |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             | Optional    |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/smrnaseq --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/smrnaseq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Optional parameters

If `--save_intermediates` is specified, the intermediate files generated in the pipeline will be saved in the output directory.

If `--save_reference` is specified, the reference files and indices generated by the pipeline will be saved to `genome/` in the output directory. When combined with `--save_intermediates`, these files are saved to `genome/` only, instead of `bowtie_index/`.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/smrnaseq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/smrnaseq releases page](https://github.com/nf-core/smrnaseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Stand-alone scripts

The `bin` directory contains some scripts used by the pipeline which may also be run manually:

- `edgeR_miRBase.r`: R script using for processing reads counts of mature miRNAs and miRNA precursors (hairpins).
  This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such a profile (such as uploading it as supplementary material for academic publications), make sure not to include cluster-specific paths to files, nor institution-specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not expose every command-line option of every tool as a pipeline parameter. For supported nf-core modules, you can pass extra tool options with the `ext.args` process directive in a Nextflow config file. This is useful, for example, when you want to change Bowtie's mismatch settings or add an option to a SAMtools command. These are module arguments, so put them in a custom config passed with `-c`; pipeline parameters belong in the CLI or a `-params-file` (see [running the pipeline](#running-the-pipeline)).

For this pipeline, the alignment module is named `BOWTIE_ALIGN`. It runs both Bowtie and `samtools view`: set `ext.args` for Bowtie options and `ext.args2` for `samtools view` options. For the separate sorting step, use `ext.args` on `SAMTOOLS_SORT`.

Genome mapping runs this module as `BOWTIE_MAP_GENOME`, and the pipeline configures it with a fully qualified selector, so a shorter selector such as `BOWTIE_ALIGN` has no effect on that step. Setting `ext.args` replaces the pipeline's defaults, so repeat any you want to keep.

For example, create `tool_args.config`:

```nextflow
process {
    withName: 'NFCORE_SMRNASEQ:GENOME_QUANT:BOWTIE_MAP_GENOME' {
        ext.args  = '-t -k 50 --best --strata -e 99999 --chunkmbs 2048 -n 0' // Bowtie: pipeline defaults, plus zero mismatches in the seed
        ext.args2 = '--no-PG'    // samtools view: example additional option
    }

    withName: SAMTOOLS_SORT {
        ext.args = '--no-PG'     // samtools sort: example additional option
    }
}
```

Then include it in your usual run command:

```bash
nextflow run nf-core/smrnaseq \
    -profile docker \
    -c tool_args.config \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh37
```

The Bowtie example sets its seed mismatch limit to zero; choose options that match your analysis and the tool's own documentation. `ext.args` values are inserted into the module command, so check the module's command and tool help before adding options that duplicate settings the pipeline already supplies. The applicable process names and argument slots can be checked in the module code; for modules running multiple tools, `ext.args`, `ext.args2`, etc. correspond to the tools in command order. See nf-core's documentation on [ext arguments](https://nf-co.re/docs/contributing/components/ext_args) and [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) for more details.

### nf-core/configs

In most cases, you will only need to create a custom config once. However, if multiple users within an organisation intend to run nf-core pipelines regularly under the same settings, then consider uploading your custom config file to the `nf-core/configs` git repository. Before uploading, ensure that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amended [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
