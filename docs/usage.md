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

`mirtrace_species` is required when using miRBase (the default), unless `mirna_species` is set (see below). Besides selecting the miRTrace species, it is used to select the miRNAs of your species from the `mature` and `hairpin` FASTA files and, unless `mirna_gtf` is given, to download the miRBase GFF3 file for mirtop (see below). Use the miRBase species code, e.g. `hsa`, `mmu` or `oar` for sheep. It cannot be set to a placeholder such as `none` to skip miRTrace, because the miRNA quantification would then find no miRNAs.

### Quantify against several species or all species (`mirna_species`, `mirgenedb_species`)

Use `--mirna_species` to filter the miRBase references to other species than `--mirtrace_species`, for example to annotate a non-model organism against several related species. It takes one species code, a comma-separated list of codes without spaces (e.g. `--mirna_species 'dre,fru,tni'`) or `all` to keep the references of all miRBase species. When it is not set, the references are filtered by `--mirtrace_species` as before. miRTrace still uses `--mirtrace_species`, which takes a single species. If `--mirtrace_species` is not set, miRTrace QC is skipped.

With a single species, the miRBase GFF3 file of that species is downloaded for mirtop unless `--mirna_gtf` is given, also when it differs from `--mirtrace_species`. With more than one species or `all`, provide a GFF3 file with the miRNA coordinates of these species via `--mirna_gtf`, since miRBase provides one GFF3 file per species, or skip mirtop with `--skip_mirtop` (see [Species without a miRBase GFF3 file](#species-without-a-mirbase-gff3-file)). The miRBase GFF3 files of the species can be concatenated into one file:

```bash
for species in dre fru tni; do
    curl -fsSL "https://mirbase.org/download/${species}.gff3"
done > mirbase_dre_fru_tni.gff3
```

```bash
nextflow run nf-core/smrnaseq \
    -profile <docker/singularity/.../institute>,<protocol> \
    --input samplesheet.csv \
    --fasta <genome.fa> \
    --mirtrace_species dre \
    --mirna_species 'dre,fru,tni' \
    --mirna_gtf mirbase_dre_fru_tni.gff3 \
    --outdir <OUTDIR>
```

With `--mirgenedb`, `--mirgenedb_species` likewise takes one species, a comma-separated list (e.g. `Hsa,Mmu`, with MirGeneDB's capitalised codes) or `all`; `--mirna_species` is not used. `--mirgenedb_gff` must then contain the coordinates of all these species, for example the MirGeneDB GFF files of the species concatenated into one file or the GFF file for all species from [MirGeneDB](https://mirgenedb.org/download). The mature and precursor FASTA files can likewise be the files for all species, as they are filtered by `--mirgenedb_species`:

```bash
cat Hsa.gff Mmu.gff > mirgenedb_hsa_mmu.gff
nextflow run nf-core/smrnaseq \
    -profile <docker/singularity/.../institute>,<protocol> \
    --input samplesheet.csv \
    --mirgenedb \
    --mirgenedb_species 'Hsa,Mmu' \
    --mirgenedb_mature mirgenedb_all_mature.fa \
    --mirgenedb_hairpin mirgenedb_all_precursor.fa \
    --mirgenedb_gff mirgenedb_hsa_mmu.gff \
    --outdir <OUTDIR>
```

When concatenating GFF files, the header lines of each file can be kept: `mirtop` reads the database name from the header of the first file and skips the other comment lines.

Please note:

- miRBase only provides GFF3 files for species with genome coordinates. Hairpins of species without entries in the GFF file are still used for the Bowtie alignments and the edgeR counts, but `mirtop` does not annotate reads on them.
- Conserved miRNAs often have identical sequences in several species. Reads of such miRNAs align equally well to the entry of each species and are reported for each of them, so the mature and hairpin counts (edgeR) and the `mirtop` isomiR tables (`joined_samples_mirtop.tsv`, `export/`) contain the same reads under several names (e.g. `dre-miR-1` and `fru-miR-1`). The `mirtop/mirna.tsv` table counts each read sequence only once, under one of these names. Sum or pick these entries in downstream analyses as needed.
- With `all`, this affects most reads of conserved miRNAs: in miRBase 22, the let-7a sequence is identical in 65 mature entries and contained in 173 hairpins. The miRNA alignments report at most 50 hits per read (`-k 50`), so for such reads only an arbitrary subset of the matching entries is counted. Prefer a list of related species over `all` where possible. The references themselves stay small (miRBase 22: about 38,600 hairpins and 48,900 mature sequences, a few MB), and building their Bowtie indices takes well under a minute. On the pipeline's test data, `all` made the hairpin alignments and `mirtop` about 2-3 times slower than a single species, raised the peak memory of `mirtop gff` to about 450 MB and of the edgeR step to about 1 GB, because of the many more multi-mapping hits and the larger count tables. The default resources cover this, but expect longer run times on full-size data.
- miRDeep2 receives the filtered mature and hairpin sequences of all selected species as known miRNAs.

### miRNA related files

Different parameters can be set for the two supported databases. By default `miRBase` will be used with the parameters below.

- `mirna_gtf`: If not supplied by the user, then `mirna_gtf` will point to the latest GFF3 file in miRBase: `https://mirbase.org/download/${params.mirtrace_species}.gff3` (the species given with `--mirna_species` if it is set to a single species). It is required when `--mirna_species` selects more than one species or `all`, unless mirtop is skipped with `--skip_mirtop`.
- `mature`: points to the FASTA file of mature miRNA sequences. Default: `https://mirbase.org/download/mature.fa`
- `hairpin`: points to the FASTA file of precursor miRNA sequences. Default: `https://mirbase.org/download/hairpin.fa`

#### Species without a miRBase GFF3 file

miRBase only provides GFF3 files for species whose miRNAs are mapped to a genome assembly. For many species in miRBase, such as sheep (`oar`) or rabbit (`ocu`), the file `https://mirbase.org/download/<species>.gff3` (for `--mirna_species`, or `--mirtrace_species` if it is not set) does not exist and the pipeline stops with an error that the miRBase GFF3 file could not be found. These species can still be analysed: their miRNAs are in the miRBase `mature` and `hairpin` FASTA files, and miRTrace supports them (`mirtrace --list-species` lists the supported species codes).

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

- `mirgenedb_gff`: The GFF file cannot be downloaded automatically due to the presence of short-term tokens in the URLs. Therefore, the user must manually provide the GFF file, either for their species of interest or for all species, by downloading it from [MirGeneDB](https://mirgenedb.org/download). The provided dataset will be automatically filtered based on the species specified with the `--mirgenedb_species` parameter, which can also be a comma-separated list of species or `all` (see [above](#quantify-against-several-species-or-all-species-mirna_species-mirgenedb_species)).
- `mirgenedb_mature`: This parameter should point to the FASTA file containing mature miRNA sequences. The file can be manually downloaded from [MirGeneDB](https://mirgenedb.org/download).
- `mirgenedb_hairpin`: This parameter should point to the FASTA file containing precursor miRNA sequences. Note that MirGeneDB does not offer a dedicated hairpin file, but the precursor sequences can be downloaded from [MirGeneDB](https://mirgenedb.org/download) and used instead.
- `mirgenedb_star` (optional): This parameter can point to the FASTA file containing the star miRNA sequences ('Star sequences' on [MirGeneDB](https://mirgenedb.org/download)). The MirGeneDB mature FASTA only contains the dominant arm of each miRNA. With this file, the star arms (named with a trailing `*`, e.g. `Hsa-Mir-3613_3p*`) are also counted in the edgeR mature tables, under the same names as in the mirtop outputs. See [miRNA identifiers](https://nf-co.re/smrnaseq/output#mirna-identifiers).

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
  - `-k` must be at least 2, because the defaults include `--strata`, which Bowtie only accepts with `-k` above 1. To report a single alignment per read, use `-M 1` instead.
  - `-M <n>` (at least 1) reports a read that has more than `n` alignments in the best stratum once, at a random one of them, with a MAPQ of 0. Reads with up to `n` alignments keep all of them (up to the `-k` limit), so `-M 1` gives exactly one alignment per aligned read. `-M` cannot be combined with `-m`, which drops such reads instead.
  - `bowtie_align_args` only changes genome mapping. To change Bowtie options, such as the number of allowed mismatches, for the mature miRNA, hairpin or seqcluster alignments, use `ext.args` in a custom config instead (see [Custom Tool Arguments](#custom-tool-arguments)).

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

### Other small RNAs: genome-wide quantification

The pipeline quantifies miRNAs only. It does not count genome-aligned reads against a genome annotation (GTF) or call small RNA clusters (see [#409](https://github.com/nf-core/smrnaseq/issues/409)). You can use the genome alignment as a starting point for this outside the pipeline. Run with `--save_intermediates` so that the BAM files are saved to `genome_quant/bam/` (see [the output documentation](https://nf-co.re/smrnaseq/output#bowtie)). Keep the following in mind:

- Only reads that did not align to the mature or hairpin miRNA references are aligned to the genome. These are mostly non-miRNA reads, such as tRNA, rRNA, snoRNA, piRNA, repeat-derived and unannotated small RNAs. With `--filter_contamination`, the classes you filter out are removed before this step, so do not filter a class you want to quantify.
- The reads are not collapsed: each alignment record comes from one sequenced read (after UMI deduplication, if used).
- Bowtie reports up to 50 alignments per read (`-k 50 --best --strata`). Bowtie 1 does not set the secondary alignment flag or the `NH` tag, and gives every alignment a MAPQ of 255. featureCounts and HTSeq use the `NH` tag to detect multi-mapping reads, so they count each alignment of a multi-mapping read as a separate, uniquely mapped read. Counting these BAM files directly against a GTF therefore inflates the counts of multi-copy classes, such as tRNAs, rRNAs and repeat-derived piRNAs.

To count reads against a GTF with [featureCounts](https://subread.sourceforge.net/featureCounts.html), use one of these approaches. The first two also change the genome alignment statistics in the MultiQC report.

- **Uniquely aligned reads only:** run the pipeline with `--bowtie_align_args '-m 1'`. Bowtie then drops reads that have more than one alignment in the best stratum.
- **One random best alignment per read:** run the pipeline with `--bowtie_align_args '-M 1'`. Bowtie gives the randomly placed alignments of multi-mapping reads a MAPQ of 0, so you can still count uniquely aligned reads only with featureCounts `-Q 1`.
- **All alignments, counted fractionally:** keep the default settings and add an `NH` tag before counting. Bowtie 1 writes the number of alignments it reported for a read in the `XM` tag, so you can rename it, for example with `samtools view -h sample.sorted.bam | perl -pe 's/\tXM:i:/\tNH:i:/' | samtools view -b -o sample.nh.bam -`, and then count with featureCounts `-M --fraction`. A read with exactly 50 alignments may have more that Bowtie did not report.

For example, to summarise the reads by gene biotype after a run with `--bowtie_align_args '-m 1'` or `'-M 1'`, where each read has one alignment:

```bash
featureCounts -s 1 -t exon -g gene_biotype -a annotation.gtf -o biotype_counts.tsv results/genome_quant/bam/*.sorted.bam
```

Do not use this command on the BAM files of a run with the default settings: it would count each of the up to 50 alignments of a read as a separate read. For a default run, add the `NH` tag as described above and count the alignments fractionally:

```bash
featureCounts -s 1 -M --fraction -t exon -g gene_biotype -a annotation.gtf -o biotype_counts.tsv *.nh.bam
```

Use `-g gene_type` for GENCODE annotations, or `-g gene_id` for counts per gene. `-s 1` counts reads on the annotated strand only, which suits most small RNA library protocols.

To find small RNA loci without relying on an annotation, use a cluster caller that handles multi-mapping reads itself. Examples are [seqcluster](https://seqcluster.readthedocs.io/) (`seqcluster prepare` and `seqcluster cluster`), which can start from the collapsed reads that the pipeline saves to `mirna_quant/seqcluster/` with `--save_intermediates`, and [ShortStack](https://github.com/MikeAxtell/ShortStack), which can start from the trimmed reads that the pipeline saves to `fastp/fastq/` with `--save_merged`.

### Contamination filtering

This step has, until now, only been tested for human data. Unexpected behaviour can occur when using it with a different species.

Contamination filtering of the sequencing reads is optional and can be invoked using the `filter_contamination` parameter. FASTA files with

- `rrna`: Used to supply a FASTA file containing rRNA contamination sequence.
- `trna`: Used to supply a FASTA file containing tRNA contamination sequence. e.g. `http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa`
- `cdna`: Used to supply a FASTA file containing cDNA contamination sequence. e.g. `ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz` The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `ncrna`: Used to supply a FASTA file containing ncRNA contamination sequence. e.g. `ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz` The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `pirna`: Used to supply a FASTA file containing piRNA contamination sequence, e.g. from [piRBase](http://bigdata.ibp.ac.cn/piRBase/) or [piRNAdb](https://www.pirnadb.org/). The FASTA file is first compared to the available miRNA sequences and overlaps are removed.
- `other_contamination`: Used to supply an additional filtering set. The FASTA file is first compared to the available miRNA sequences and overlaps are removed.

How the filter works:

- The adapter-trimmed reads are aligned with Bowtie2 (`--very-sensitive-local -k 1`) to each database in turn, in the fixed order rRNA, tRNA, cDNA, ncRNA, piRNA, other. Only the reads that do not align are passed on to the next database, so each read is counted for the first class it aligns to. The order matters when databases overlap (for example, rRNA and tRNA sequences are also part of the Ensembl ncRNA FASTA).
- Before aligning, the cDNA, ncRNA, piRNA and other databases are compared to the miRNA hairpins with BLAT, and every entry with a significant hit (e-value < 1e-5) is removed from the database, so that miRNA reads are not filtered out. The rRNA and tRNA databases are used as given.
- The reads that align to a database are removed from the analysis, and only their number is kept: their alignments and sequences are not saved to the output directory. The remaining reads are used for miRNA quantification, genome mapping and miRDeep2.
- The number of reads removed per class is shown in the "Contamination Filtering" plot of the MultiQC report. With `--save_intermediates`, the per-sample numbers (`*.contamination_mqc.yaml`) and the filtered reads (`*.filtered.fastq.gz`) are also saved to `contaminant_filter/filter/` (see [the output documentation](https://nf-co.re/smrnaseq/output#bowtie2)).

### Small RNA biotypes other than miRNAs

nf-core/smrnaseq quantifies miRNAs: it produces count tables per mature miRNA, per hairpin and per isomiR. It does not produce count tables per tRNA, tRNA fragment, piRNA, rRNA, snoRNA or other small RNA. Data for these biotypes is available at the following levels:

| What                                                                           | Where                                                                                                                                      | Parameters                                                       |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------- |
| Fraction of reads that are miRNA, rRNA, tRNA, artifacts or unknown (QC)        | miRTrace "RNA type" plot in MultiQC, `mirtrace/<sample>/`                                                                                  | `--mirtrace_species` (species supported by miRTrace)             |
| Number of reads per database (rRNA, tRNA, cDNA, ncRNA, piRNA, other)           | "Contamination Filtering" plot in MultiQC, `contaminant_filter/filter/*.contamination_mqc.yaml`                                            | `--filter_contamination` plus one or more of the databases above |
| Genome alignments of all reads that did not align to a mature miRNA or hairpin | `genome_quant/bam/*.sorted.bam` (BAM files only with `--save_intermediates`); samtools `stats`, `flagstat` and `idxstats` are always saved | `--fasta` (and optionally `--bowtie_index`)                      |

The contamination filter reports one total per database, not counts per sequence. The `idxstats` files of the genome mapping count reads per chromosome, not per gene.

To get counts per small RNA (for example per tRNA gene or per piRNA locus), count the genome alignments against an annotation of your choice, outside the pipeline. [Other small RNAs: genome-wide quantification](#other-small-rnas-genome-wide-quantification) describes how to do this with featureCounts, including how to handle multi-mapping reads. In addition:

- Do not set `--filter_contamination`, or at least do not set the databases of the biotypes you want to count: reads removed by the filter are not mapped to the genome.
- Use an annotation that matches your genome assembly, for example the [GENCODE](https://www.gencodegenes.org/) or [Ensembl](https://www.ensembl.org/) GTF (snoRNA, snRNA, rRNA, misc_RNA and more, see the `gene_type`/`gene_biotype` attribute), [GtRNAdb](http://gtrnadb.ucsc.edu/) (tRNA genes), [piRBase](http://bigdata.ibp.ac.cn/piRBase/) (piRNAs) or [RNAcentral](https://rnacentral.org/) (genome coordinates of many ncRNA databases). featureCounts reads GTF/GFF or its own SAF format, so convert BED annotations (as offered by GtRNAdb) to SAF first.
- Reads that align to a mature miRNA or hairpin are not part of the genome BAM files: take miRNA counts from `mirna_quant/`.
- Mature tRNAs carry a non-templated 3' `CCA` and lack their introns, and tRNA-derived fragments are short. Genome alignment misses many of these reads. Dedicated tools that align against mature tRNA sequences give better tRNA and tRNA fragment counts.

If you need full annotation and quantification of all small RNA classes, consider using a dedicated tool on the adapter-trimmed reads (saved to `fastp/fastq/` with `--save_merged`), such as [unitas](https://www.smallrnagroup.uni-mainz.de/software.html), [SPORTS](https://github.com/junchaoshi/sports1.1), [exceRpt](https://github.com/gersteinlab/exceRpt) or, for tRNA fragments, [MINTmap](https://github.com/TJU-CMC-Org/MINTmap). Support for further biotypes in the pipeline is discussed in [#548](https://github.com/nf-core/smrnaseq/issues/548).

## mirDeep2

[miRDeep2](https://github.com/rajewsky-lab/mirdeep2) predicts novel miRNAs from the reads and also reports how well the data supports the known miRNAs. It is the only step of the pipeline that looks for miRNAs that are not in the miRBase or MirGeneDB reference. The pipeline does not detect or quantify other small RNA classes (for example piRNAs) on their own; the `pirna` FASTA is only used to remove reads in [contamination filtering](#contamination-filtering).

### Requirements and inputs

- miRDeep2 runs by default whenever the pipeline has a genome FASTA, and uses the genome Bowtie 1 index to map the reads. Without a genome FASTA, the pipeline skips genome quantification and miRDeep2 and prints a warning at start-up. Use `--skip_mirdeep` to turn miRDeep2 off when a genome is given. See [Genome quantification and miRDeep2](#genome-quantification-and-mirdeep2) for how to give the genome and the index.
- Reads: the trimmed reads of each sample, after UMI handling and, with `--filter_contamination`, after contamination filtering. All of these reads go to miRDeep2, including reads that map to known miRNAs. The pipeline converts them to FASTA, then collapses and maps them with `mapper.pl -c -j -m -v`.
- Known miRNAs: the mature and hairpin sequences of your species, taken from `--mature` and `--hairpin` (or `--mirgenedb_mature`, `--mirgenedb_hairpin` and, if given, `--mirgenedb_star`) and filtered by `--mirtrace_species` (or `--mirgenedb_species`). Mature miRNAs of related species are not passed to miRDeep2, so novel candidates get no seed conservation support from other species.
- miRDeep2 runs on each sample separately, single-threaded (`process_single`). By default, results are not combined across samples. With `--mirdeep2_consensus_min_samples`, the pipeline also writes the loci found in several samples (see [Cross-sample consensus of miRDeep2 loci](#cross-sample-consensus-of-mirdeep2-loci)). Novel candidates are not added to the miRNA count tables (`mirna_quant/`). To count novel miRNAs in all samples, add the consensus precursor and mature sequences of the candidates you trust to copies of your `--hairpin` and `--mature` FASTA files and run the pipeline again. Give them names with your species prefix, for example `hsa-novel-1`, because the pipeline keeps only entries that match `<species>-`. They are then counted in the mature and hairpin count tables (`mirna_quant/edger_qc/`). mirtop only annotates miRNAs that are in the `--mirna_gtf` file.

See [the output documentation](https://nf-co.re/smrnaseq/output#mirdeep2) for how to read the results.

### Tuning miRDeep2

By default, the pipeline passes these options to `miRDeep2.pl`:

- `-P` (only with miRBase references, not with `--mirgenedb`): current miRBase releases (v18 and later) name mature miRNAs with `-5p`/`-3p` suffixes, and the miRDeep2 documentation asks for `-P` with such names. MirGeneDB names end in `_5p`/`_3p`, which `-P` does not recognise, so the pipeline leaves it out in MirGeneDB mode. `-P` only changes how miRDeep2 splits known miRNAs into 5p and 3p arms in its internal quantification report, which the pipeline does not publish. The `result_*` files are the same with or without it.
- `-d`: skips the per-miRNA structure PDFs. miRDeep2 draws one PDF of about 0.3 MB for every reported miRNA and every known miRNA with reads, which can add up to hundreds of MB per sample. Because of `-d`, the miRNA IDs in the HTML report link to PDFs that do not exist.

You can change the options with `ext.args` (see [Custom Tool Arguments](#custom-tool-arguments)). Setting `ext.args` replaces the defaults, so keep `-P` when you use miRBase references. For example, to get the structure PDFs, leave out `-d`:

```groovy
process {
    withName: 'MIRDEEP2_MIRDEEP2' {
        ext.args = '-P'
    }
}
```

The pipeline then publishes the PDFs to `mirdeep2/pdfs_<sample>/`, and the links in the HTML report work. Other useful options:

- `-b <score>`: only report novel candidates with at least this miRDeep2 score (default 0).

### Exit status 255

If the software encounters an error with exit status 255, it will be ignored, and the pipeline will continue to complete. In such cases, the pipeline will log a note that includes the path to the work directory where the issue occurred. You can inspect this work directory to examine your input data and troubleshoot the issue.

Error 255 is typically related to the core algorithm of miRDeep generating empty output files. This often happens when the reads being processed do not correspond to putative mature miRNA sequences, or if the provided precursors do not meet the criteria for valid miRNA precursors, both of which may stem from the input reads used. A common cause of this error is running the pipeline with a small subset of the input reads.

### Cross-sample consensus of miRDeep2 loci

miRDeep2 reports known and novel miRNAs for each sample separately, and the provisional IDs of novel miRNAs differ between samples. To get one set of loci across samples, use `--mirdeep2_consensus_min_samples <N>`: loci of the per-sample BED files that overlap by at least one base on the same strand are merged, and loci reported by at least `N` samples are kept (for example `1` for the union, or the number of samples for loci found in all of them). Loci that only touch, without sharing a base, are not merged. The result is written to `mirdeep2/mirdeep2_consensus.bed`, a standard six-column BED file for genome browsers, and to `mirdeep2/mirdeep2_consensus.tsv`, which lists the samples, the highest miRDeep2 score and probability and the per-sample miRDeep2 IDs of each locus (see [the output documentation](https://nf-co.re/smrnaseq/output#cross-sample-consensus-of-mirdeep2-loci)). The default `0` disables this step.

Add `--mirdeep2_consensus_min_probability <P>` to use only miRNAs whose estimated true-positive probability in the per-sample `result_{sample.id}.csv` is at least `P` percent (for example `75`; `91 +/- 2%` counts as 91). Entries below the threshold are removed in each sample before the samples are merged, so they do not count towards `--mirdeep2_consensus_min_samples`. miRDeep2 gives no probability (`-`) for known miRBase miRNAs in a sample in which it found no novel miRNA candidates. Such entries are removed when the filter is on, as they have no estimate that reaches the threshold. The default `0` disables the filter and keeps all entries. The per-sample miRDeep2 output files are not changed.

### UMI handling

With `--with_umi`, the pipeline removes PCR duplicates using unique molecular identifiers (UMIs) before any alignment. The steps run in this order:

1. `umi_tools extract` runs on the raw reads. It uses `--umitools_extract_method` and `--umitools_bc_pattern` to move the UMI into the read name and to remove the sequences around it, so that only the insert of the RNA molecule remains. Reads that do not match the pattern are discarded and counted in the `umi_tools extract` log.
2. `fastp` filters the inserts by quality and length (`--fastp_min_length`, `--fastp_max_length`) and clips `--clip_r1` bases from the 5' end; `--three_prime_clip_r1` bases are then clipped from the 3' end. With `--umitools_extract_method regex`, the pattern is expected to remove the kit adapter, so `fastp` does not trim adapters. With `string`, `fastp` trims `--three_prime_adapter` as without UMIs. The pipeline's list of known miRNA adapters is not used with UMIs.
3. The reads are deduplicated on **UMI + insert sequence** with [UMICollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse) in FastQ mode. UMICollapse cannot read the UMI from the read name, so the pipeline moves the UMI to the start of the read sequence, runs UMICollapse, and removes the UMI from the sequence again. By default, reads whose UMI + insert sequences differ by one substitution are merged (see [Choosing a UMI strategy](#choosing-a-umi-strategy)). Reads of different lengths are never merged.

The deduplicated inserts are then used for all later steps; the UMI stays in the read name. miRTrace QC runs on the raw reads as without UMIs.

Step 3 can be replaced by an experimental alignment-based deduplication with `--umi_dedup_method fgumi` (see [Experimental: alignment-based deduplication with fgumi](#experimental-alignment-based-deduplication-with-fgumi)). The default, `--umi_dedup_method umicollapse`, is the workflow described above.

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
- `--umitools_method`: despite its name, this sets the UMICollapse deduplication algorithm (`--algo`): `dir` (directional, default), `adj` (adjacency) or `cc` (connected components). All three merge UMIs that differ by one substitution. See [UMI deduplication algorithm](#umi-deduplication-algorithm). It is not used with `--umi_dedup_method fgumi`.
- `--umi_dedup_method`: `umicollapse` (default) or `fgumi` (experimental). See [Experimental: alignment-based deduplication with fgumi](#experimental-alignment-based-deduplication-with-fgumi).
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

#### Memory use of UMI deduplication

In `fastq` mode, UMICollapse reads the whole FASTQ file of a sample into memory before it collapses reads. As the UMI makes most reads distinct, memory use and run time grow with the sequencing depth of a sample. The `UMICOLLAPSE_FASTQ` process requests 200 GB of memory by default (within your [maximum resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources)), and the Java heap is set to 90% of the memory of the task. A Java `OutOfMemoryError` is not retried automatically. For very deep samples, give this process more memory and time in a custom config passed with `-c`:

```groovy
process {
    withName: '.*:UMICOLLAPSE_FASTQ' {
        memory = 300.GB
        time   = 48.h
    }
}
```

#### Sequence-based and alignment-based UMI deduplication

Many UMI workflows for mRNA deduplicate after alignment instead: reads are grouped by alignment position (reference, strand and 5' start), and only the UMIs of reads at the same position are compared (e.g. `umi_tools dedup`, UMICollapse BAM mode or `fgumi dedup`). The pipeline offers this as an experimental option, `--umi_dedup_method fgumi` (see [Experimental: alignment-based deduplication with fgumi](#experimental-alignment-based-deduplication-with-fgumi)). Neither approach is wrong for small RNA, but they define a duplicate differently:

- **Sequence-based (pipeline default)**: needs no reference and keeps isomiRs apart. Because each insert sequence and length is its own group, abundant miRNAs are split over several groups, which reduces UMI saturation (see [Choosing a UMI strategy](#choosing-a-umi-strategy)). The one allowed substitution (default `-k 1`) covers the UMI and the insert together: a read with a single sequencing error, in the UMI or in the insert, is still collapsed with its duplicate. Only reads with two or more errors (for example one in the UMI and one in the insert) are kept as separate molecules, and reads of different length are never compared.
- **Alignment-based**: sequencing errors in the insert do not matter, and memory use is lower because only the reads at one position are compared. However, almost all reads of a miRNA start at the same position, so all its 3' isomiRs are pooled and told apart by their UMI alone. For highly expressed miRNAs this increases UMI saturation and can compress their counts when UMIs that differ by one base are merged. The pipeline's quantification alignments use `-k 50`, so deduplicating them directly would treat reads that align to several miRNAs or loci separately at each position, and the same read could be kept at one locus and removed at another. `--umi_dedup_method fgumi` therefore uses its own alignment with one record per read.

To deduplicate the pipeline's own alignments by position instead, you can extract the UMIs before running the pipeline (with `umi_tools extract` and the same `--bc-pattern`, which moves the UMI to the read name and removes the kit sequences around it), run the pipeline without `--with_umi` and with `--save_intermediates`, and then deduplicate the `*.sorted.bam` files in `mirna_quant/bam/` and `genome_quant/bam/` with `umi_tools dedup` or `umicollapse bam`, followed by `samtools idxstats` for counts. Note that the mirtop, edgeR and miRDeep2 results of such a run are not deduplicated.

#### Experimental: alignment-based deduplication with fgumi

> [!WARNING]
> `--umi_dedup_method fgumi` is experimental. It is opt-in, the pipeline prints a warning at start-up when it is used, and its behaviour, defaults and output files may change in a future release. The default UMI workflow (`--umi_dedup_method umicollapse`) is not affected.

With `--with_umi --umi_dedup_method fgumi`, steps 1 and 2 of the UMI workflow stay the same, and step 3 (UMICollapse) is replaced by a deduplication by alignment position with [fgumi](https://github.com/fulcrumgenomics/fgumi):

1. The extracted and filtered inserts are aligned with Bowtie 1 (`-k 1 --best`, plus `--norc` with `--mirna_forward_only`) to a dedicated deduplication reference (see below). Each read gets at most one alignment.
2. `fgumi copy-umi` copies the UMI from the read name into the `RX` tag, and `fgumi sort` sorts the alignments in template-coordinate order.
3. `fgumi dedup --strategy identity --remove-duplicates --include-unmapped` keeps one read (the one with the highest base qualities) for each group of reads with the same reference sequence, 5' position, strand and **identical UMI**.
4. `fgumi fastq` writes the remaining reads back to FASTQ, in their original orientation and with the UMI still in the read name. These reads are used for all later steps (contamination filter, mature and hairpin quantification, mirtop, genome quantification, miRDeep2), exactly like the UMICollapse output.

```bash
nextflow run nf-core/smrnaseq \
    -profile docker,qiaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --mirtrace_species hsa \
    --with_umi \
    --umi_dedup_method fgumi
```

**Deduplication reference ("windows").** The reference is built from the mature and hairpin FASTA files of the selected species (miRBase or MirGeneDB, after the same species filter as for quantification, including `--mirgenedb_star`). It has one record per distinct mature miRNA sequence: the mature sequence plus 5 nt on each side, taken from the first hairpin in the hairpin FASTA that contains it. The flanks are shortened at the ends of the hairpin, and a mature sequence that is found in no hairpin is used without flanks. Mature miRNAs with the same sequence (for example hsa-let-7a-5p, which comes from let-7a-1, let-7a-2 and let-7a-3) share one record. The reference and a table of its records are saved to `umi_dedup/reference/`. This design follows the benchmark in [#535](https://github.com/nf-core/smrnaseq/issues/535#issuecomment-7086624082):

- With the hairpins as reference, identical reads are distributed at random over paralogous hairpins and escape deduplication, which inflated the counts of multi-copy miRNAs by 10–60% on real data.
- With the mature sequences as reference, isomiRs that are longer than the mature sequence do not align and are not deduplicated (15–21% count error in simulation).

The number of flanking nucleotides can be changed with `ext.args = '-v flank=<n>'` for the `UMI_DEDUP_WINDOWS` process.

**What you should know before using it:**

- **Reads that do not align to the windows reference are not deduplicated.** They are passed through unchanged, so that they are still available for the later steps. This affects everything that is not a known miRNA of the selected species: other small RNAs, novel miRNAs and reads with many mismatches. The hairpin, genome quantification and miRDeep2 results are therefore only partly deduplicated. `passthrough_templates` in the metrics file reports how many reads this concerns.
- **IsomiRs with the same 5' start and the same UMI are pooled.** Reads are grouped by position and strand, not by sequence: two different molecules of one miRNA (for example two 3' isomiRs) that carry the same UMI by chance are counted once. This matters for the most abundant miRNAs, where a large part of all possible UMIs is in use.
- **UMIs are matched exactly** (`--strategy identity`). A read with a PCR or sequencing error in its UMI is counted as a separate molecule: a near-uniform overcount of about 5–9% in the benchmark. In return, distinct molecules with similar UMIs are not merged at highly expressed miRNAs (see [Choosing a UMI strategy](#choosing-a-umi-strategy)). `--umitools_method` only sets the UMICollapse algorithm and has no effect here, and no UMI saturation statistics or saturation warnings are produced.
- **Reads whose UMI contains an `N` are removed** by `fgumi dedup` (`filtered_ns_in_umi` in the metrics file).
- **Accuracy.** In the benchmark, the per-miRNA count error was 0.6–1.1% for the default UMICollapse workflow and 1.3–2.5% for fgumi with a windows reference. The benchmark used a prototype of the reference; the figures were not measured again with the implementation in the pipeline.
- **Cost.** The method adds one Bowtie 1 alignment of all reads to a small reference and four fgumi steps for each sample. `fgumi dedup` itself is fast and its memory use is bounded (5–23 s CPU for 4.5M reads in the benchmark, compared with about 440 s for UMICollapse `dir`), so the method can help for deep samples where UMICollapse needs a lot of memory ([#530](https://github.com/nf-core/smrnaseq/issues/530)). Most of the remaining time is `umi_tools extract`.
- **`umi_tools extract` is still used** to extract the UMI. `fgumi extract` describes reads with positional read structures and cannot search for a kit sequence, which the QIAseq read structure needs (the UMI follows a variable-length insert and the kit adapter).
- MultiQC does not show the fgumi metrics yet. They are saved to `umi_dedup/fgumi/` (see the [output documentation](https://nf-co.re/smrnaseq/output#fgumi-deduplicate-experimental)).

The fgumi and Bowtie options are set with `ext.args` in `conf/modules.config` (processes `BOWTIE_MAP_UMI_DEDUP`, `FGUMI_COPYUMI`, `FGUMI_SORT`, `FGUMI_DEDUP`) and can be changed in a custom config (see [Custom Tool Arguments](#custom-tool-arguments)). Keep `-k 1` for Bowtie: `fgumi dedup` needs one record per read.

### Choosing a UMI strategy

UMI deduplication has to decide when two reads come from the same molecule. Reads with the same insert and UMI are duplicates. Reads whose UMIs differ by one base can be duplicates with a PCR or sequencing error in the UMI, or two different molecules whose UMIs happen to be similar. The default UMI workflow (`--umi_dedup_method umicollapse`) offers two strategies:

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

**Alignment-based alternative (experimental).** `--umi_dedup_method fgumi` deduplicates by alignment position with [fgumi](https://github.com/fulcrumgenomics/fgumi) and exact UMI matching (see [Experimental: alignment-based deduplication with fgumi](#experimental-alignment-based-deduplication-with-fgumi)).

**Planned.** [chelae](https://github.com/fulcrumgenomics/chelae) as an optional fast UMI extractor and the fgumi metrics in MultiQC (which needs a MultiQC release with the fgumi module) are planned for a later release ([#535](https://github.com/nf-core/smrnaseq/issues/535)). An alignment-based mode for UMICollapse (BAM mode) is proposed in [#619](https://github.com/nf-core/smrnaseq/issues/619).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 2 columns ("sample" and "fastq_1"), and a header row as shown in the examples below.

If a second fastq file is provided using another column, the extra data are ignored by this pipeline. The smRNA species should be sufficiently contained in the first read, and so the second read is superfluous data in this smRNA context.

```bash
--input '[path to samplesheet file]'
```

A samplesheet with a single sample is supported, for example for a first test run:

```console
sample,fastq_1
SAMPLE1,/path/to/SAMPLE1.fastq.gz
```

> [!NOTE]
> nf-core/smrnaseq v2.3.1 and earlier failed in miRTrace with `Cannot cast object '<file>.fastq.gz' with class 'nextflow.processor.TaskPath' to class 'java.util.List'` when the samplesheet contained only one sample ([#249](https://github.com/nf-core/smrnaseq/issues/249)). If you see this error, update to the latest release, or add a second sample to the samplesheet.

### Multiple runs of the same sample

The `sample` identifiers should match between runs of resequenced samples. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet must have at least 2 columns (`sample` and `fastq_1`). A third column can be added if the sample is paired-end (`fastq_2`).

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

| Column    | Description                                                                                                                                       | Requirement |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Sample names cannot contain spaces. | Mandatory   |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                        | Mandatory   |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                        | Optional    |

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

The pipeline runs this module in four Bowtie alignment steps. Each is configured with a fully qualified selector, so a shorter selector such as `BOWTIE_ALIGN` has no effect on them; use the selector of the step you want to change:

| Step                                                                          | Selector                                            |
| ----------------------------------------------------------------------------- | --------------------------------------------------- |
| Mature miRNA alignment                                                        | `NFCORE_SMRNASEQ:MIRNA_QUANT:BOWTIE_MAP_MATURE`     |
| Hairpin alignment of the reads that did not align to mature miRNAs            | `NFCORE_SMRNASEQ:MIRNA_QUANT:BOWTIE_MAP_HAIRPIN`    |
| Hairpin alignment of the collapsed reads, used by mirtop                      | `NFCORE_SMRNASEQ:MIRNA_QUANT:BOWTIE_MAP_SEQCLUSTER` |
| Genome alignment of the reads that did not align to mature miRNAs or hairpins | `NFCORE_SMRNASEQ:GENOME_QUANT:BOWTIE_MAP_GENOME`    |

All four steps use the same basic Bowtie defaults: `-t -k 50 --best --strata -e 99999 --chunkmbs 2048`. With [`--mirna_forward_only`](#mirna-reference-orientation), which `-profile qiaseq` sets, the mature, hairpin and seqcluster steps also get `--norc`. For genome mapping only, you can append options with the [`bowtie_align_args`](#genome) parameter instead of a custom config. Setting `ext.args` replaces the pipeline's defaults for that step, including any options added there by pipeline parameters, so repeat any you want to keep.

For example, to allow at most one mismatch in the mature miRNA and hairpin alignments and add example options for SAMtools, create `tool_args.config`:

```nextflow
process {
    withName: 'NFCORE_SMRNASEQ:MIRNA_QUANT:BOWTIE_MAP_MATURE' {
        ext.args  = '-t -k 50 --best --strata --chunkmbs 2048 -v 1' // Bowtie: pipeline defaults, but at most 1 mismatch in the whole read
        // With --mirna_forward_only or -profile qiaseq, add --norc here as well, e.g. '... -v 1 --norc'
        ext.args2 = '--no-PG'    // samtools view: example additional option
    }

    withName: 'NFCORE_SMRNASEQ:MIRNA_QUANT:BOWTIE_MAP_HAIRPIN' {
        ext.args  = '-t -k 50 --best --strata --chunkmbs 2048 -v 1' // add --norc with --mirna_forward_only or -profile qiaseq
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

This example replaces the whole `ext.args` of the two steps, so it also drops the `--norc` that `--mirna_forward_only` (set by `-profile qiaseq`) adds there: without `--norc` in your `ext.args`, reads are again aligned to both strands of the miRNA references. If you use `--mirna_forward_only` or `-profile qiaseq`, append `--norc` to the `ext.args` of every mature, hairpin and seqcluster step you override.

Bowtie 1 has two alignment modes. In the default [`-n` mode](https://bowtie-bio.sourceforge.net/manual.shtml#the--n-alignment-mode), `-n` limits the mismatches in the seed at the start of the read (default 2 in the first `-l` 28 bases), and `-e` caps the sum of the base qualities at all mismatched positions. The pipeline sets `-e 99999`, so mismatches after the seed are effectively unlimited. In [`-v` mode](https://bowtie-bio.sourceforge.net/manual.shtml#the--v-alignment-mode), Bowtie allows at most that many mismatches over the whole read and ignores `-l`, `-e` and base qualities, so the example leaves out `-e 99999`. `-v` cannot be combined with `-n`. Stricter settings in an earlier step pass more reads on to the next one (mature, then hairpin, then genome), so they can change downstream counts too. Choose options that match your analysis and the tool's own documentation. `ext.args` values are inserted into the module command, so check the module's command and tool help before adding options that duplicate settings the pipeline already supplies. The applicable process names and argument slots can be checked in the module code; for modules running multiple tools, `ext.args`, `ext.args2`, etc. correspond to the tools in command order. See nf-core's documentation on [ext arguments](https://nf-co.re/docs/contributing/components/ext_args) and [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) for more details.

### nf-core/configs

In most cases, you will only need to create a custom config once. However, if multiple users within an organisation intend to run nf-core pipelines regularly under the same settings, then consider uploading your custom config file to the `nf-core/configs` git repository. Before uploading, ensure that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amended [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Nextflow plugins

The pipeline loads these Nextflow plugins (versions are pinned in `nextflow.config`). Nextflow downloads them from the [Nextflow plugin registry](https://registry.nextflow.io) at start-up; on systems without internet access, install them beforehand, for example with `nextflow plugin install nf-co2footprint@1.4.1`.

| Plugin                                                            | What it does in this pipeline                                                                                                         |
| ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| [nf-schema](https://github.com/nextflow-io/nf-schema)             | Validates parameters and the samplesheet, and prints `--help`.                                                                        |
| [nf-core-utils](https://github.com/nf-core/nf-core-utils)         | Start-up checks, the software versions file, the methods description with citations in the MultiQC report and the completion email.   |
| [nf-prov](https://github.com/nextflow-io/nf-prov)                 | Writes a [BioCompute Object](https://biocomputeobject.org/) provenance report to `pipeline_info/`.                                    |
| [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) | Estimates the energy use and CO₂-equivalent emissions of the run and writes a report, a summary and a trace file to `pipeline_info/`. |

### CO₂ footprint estimates

[nf-co2footprint](https://nextflow-io.github.io/nf-co2footprint) estimates the energy use of every task from its resource usage and converts it to CO₂-equivalent emissions with the carbon intensity of the electricity, following the [Green Algorithms](https://www.green-algorithms.org) method. By default the pipeline uses a global average carbon intensity and makes no network requests. For estimates that match where you run the pipeline, set your [Electricity Maps zone](https://portal.electricitymaps.com/docs/getting-started#geographical-coverage) in a custom config (`-c`):

```groovy title="co2footprint.config"
co2footprint {
    location = 'DE' // your zone code
}
```

This uses the 2024 yearly average carbon intensity for that zone. For real-time values, also give an [Electricity Maps](https://www.electricitymaps.com/) API key. Store it as a [Nextflow secret](https://www.nextflow.io/docs/latest/secrets.html) (`nextflow secrets set EM_API_KEY <key>`) rather than in the config file:

```groovy title="co2footprint.config"
co2footprint {
    location = 'DE'
    emApiKey = secrets.EM_API_KEY
}
```

The plugin documentation lists further options, such as the power usage effectiveness (PUE) of your data centre and the power draw of your CPUs.

### Provenance reports

[nf-prov](https://github.com/nextflow-io/nf-prov) writes a BioCompute Object (`pipeline_info/manifest_<timestamp>.bco.json`) at the end of every run. It can also write a [Workflow Run RO-Crate](https://www.researchobject.org/workflow-run-crate/), which describes the run together with its inputs and published results. The crate is not written by default:

- it copies the pipeline's `README.md`, `main.nf` and `nextflow_schema.json`, local input files and the **fully resolved Nextflow configuration** into the crate directory. The resolved configuration contains every value from your config files, including credentials or tokens set there (also Nextflow secrets that a config file reads). Check your configuration before you enable the crate, and do not enable it if the results are shared and your configuration holds secrets.
- to describe the published results, the crate must be the output directory itself, so these files are written to the top of `--outdir`.

To write the crate, add it with a custom config:

```groovy title="ro-crate.config"
prov {
    formats {
        wrroc {
            file      = "${params.outdir}/ro-crate-metadata.json"
            overwrite = true
            license   = 'https://spdx.org/licenses/CC-BY-4.0' // the licence of your results
        }
    }
}
```

The [WRROC options](https://github.com/nextflow-io/nf-prov/blob/master/docs/WRROC.md) also let you name the person and organisation that ran the pipeline.

### Turning plugin reports off

The reports are configured in `nextflow.config` and can be switched off with a custom config:

```groovy title="no_plugin_reports.config"
prov.enabled = false // no BioCompute Object
co2footprint {
    trace.enabled      = false
    summary.enabled    = false
    report.enabled     = false
    provenance.enabled = false
}
```

The completion email is sent by nf-core-utils with Nextflow's [mail settings](https://www.nextflow.io/docs/latest/notifications.html#mail-configuration). For runs launched from Seqera Platform, which has its own notifications, the email is only sent with `--force_completion_email`.

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
