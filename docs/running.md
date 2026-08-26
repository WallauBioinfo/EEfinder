# Running EEfinder

EEfinder exposes a **single command**. It takes a genome FASTA plus the
protein database, its metadata CSV and the host-gene baits, and produces the
endogenous element sequences, their taxonomy table and their flanking regions.

## Example run (with the bundled `test_files/`)

The repository ships a small example dataset in
[`test_files/`](https://github.com/WallauBioinfo/EEfinder/tree/master/test_files)
so you can try the full pipeline end-to-end:

| File | Role | CLI flag |
|------|------|----------|
| `Ae_aeg_Aag2_ctg_1913.fasta` | Query genome contig (*Aedes aegypti* Aag2) | `-in` |
| `virus_subset.fa` | Viral protein database | `-db` |
| `virus_subset.csv` | Metadata for the viral database | `-mt` |
| `filter_subset.fa` | Host-gene bait proteins | `-bt` |

From the repository root, with the `EEfinder` environment active:

```bash
eefinder \
  -in test_files/Ae_aeg_Aag2_ctg_1913.fasta \
  -od results_test \
  -db test_files/virus_subset.fa \
  -mt test_files/virus_subset.csv \
  -bt test_files/filter_subset.fa \
  -ln 1000 \
  -id \
  -p 2 \
  -lm 100
```

- `-ln 1000` lowers the minimum contig length (the example contig is shorter than
  the 10000 nt default, which would otherwise filter it out).
- `-id` builds the BLAST indexes for the two databases (needed on the first run
  against a given database).
- `-p 2` uses two threads.
- `-lm 100` merges same-taxon elements within 100 nt.

## Running with DIAMOND

`-md`/`--mode` selects the search engine. The default is `blastx`; the other six
values are DIAMOND sensitivity settings, in increasing order of sensitivity and
runtime: `fast`, `mid-sensitive`, `sensitive`, `more-sensitive`,
`very-sensitive`, `ultra-sensitive`.

```bash
eefinder -in <genome.fasta> -db <proteins.fa> -mt <proteins.csv> \
         -bt <host_baits.fa> -od results/ -md fast -id
```

`-md` selects the engine for **both** searches and for the `-id` index build, so
a database indexed for BLAST must be re-indexed before its first DIAMOND run
(and vice versa).

(diamond-sensitivity)=
### The sensitivity trade-off

```{important}
DIAMOND is faster than BLAST, but **at a real cost in sensitivity**, and the
loss falls exactly where endogenous-element studies operate: on the most
divergent sequences.
```

DIAMOND seeds its alignments on a reduced amino-acid alphabet. Its authors
report that this does not compromise sensitivity, but the EEfinder benchmark
found a clear effect on the highly divergent sequences typical of EVE studies.
Screening the *Aedes aegypti* Aag2 genome (GCA_021653915) against a viral
protein database gave:

| Search mode | Elements recovered | Identity range | Runtime |
|-------------|--------------------|----------------|---------|
| `blastx` (BLAST) | **481** | 10.2–100 % | 56 h 14 min |
| `-md very-sensitive` (DIAMOND) | 225 | 16.2–100 % | 42 h 34 min |
| `-md fast` (DIAMOND) | 126 | 16.9–100 % | 8 h 6 min |

DIAMOND's most sensitive setting recovered **less than half** the elements BLAST
did, for a runtime saving of only ~25 %; `fast` mode recovered about a quarter of
them, but seven times faster. The identity ranges show where the difference
comes from — BLAST reaches down to 10.2 % identity, while both DIAMOND modes
floor out around 16 %, i.e. the oldest and most degraded integrations are the
ones that go missing.

**Use `blastx` (the default) whenever sensitivity matters**, which for EE
discovery is almost always. The DIAMOND modes are appropriate for rapid
exploratory screens, for very large sample sets, or when only comparatively
recent, well-conserved integrations are of interest — with the caveat that the
resulting element counts are not comparable to BLAST-based ones.

The full benchmark is reported in section 3.2 (*Benchmark of alignment tools*) of
[Dias, Dezordi & Wallau (2024)](https://doi.org/10.1016/j.csbj.2024.10.012)
([PMC11532726](https://pmc.ncbi.nlm.nih.gov/articles/PMC11532726/)).

```{note}
The DIAMOND modes also run with the subprocess stderr routed to `/dev/null`, so
a misbehaving `diamond` build can fail silently — verify the binary if a DIAMOND
run produces no hits.
```

## Pipeline steps

Each step is a small **side-effect class**; files flow through disk with
accreting suffixes (`.rn`, `.fmt`, `.blastx`, `.filtred`, `.bed`, `.tax`, …), so
every intermediate can be inspected after the run (see
[Outputs](output.md#intermediate-files-tmp_files)).

1. **InsertPrefix** — prefix every FASTA header with `>PREFIX/…` so element
   names carry the run identity (`{prefix}.rn`).
2. **RemoveShortSequences** — drop contigs below `--length` (`{prefix}.rn.fmt`).
3. **MakeDB** — build the BLAST or DIAMOND databases for `-db` and `-bt`
   (only with `--index_databases`).
4. **SimilaritySearch** — the main search of the genome against the protein
   database. `blastx` runs with `-evalue 1e-5 -word_size 3 -matrix BLOSUM45
   -max_intron_length 100 -soft_masking true`; DIAMOND with `-e 1e-5 --matrix
   BLOSUM45 -k 500 --max-hsps 0 --<mode>`.
5. **FilterTable** — collapse redundant hits per query and range (see
   [Range junction](#range-junction-arg)),
   record the sense of each hit, and drop alignments shorter than **33 aa**.
6. **GetFasta** — extract the putative EE sequences with `bedtools getfasta`.
7. **SimilaritySearch + FilterTable + CompareResults** — search the putative EEs
   against the host-gene baits and remove every element whose best host-bait hit
   scores higher than its best database hit.
8. **GetTaxonomy** — left-join the surviving hits to the metadata CSV on the
   protein accession.
9. **GetAnnotBed → MergeBed → RemoveAnnotation → GetFasta** — merge truncated
   fragments of the same taxon and sense within `--limit` nt, where "same taxon"
   is decided at `--merge_level`, then re-extract the merged sequences.
10. **MaskClean** *(optional)* — with `--clean_masked`, drop elements whose
    soft-masked/ambiguous base content exceeds `--mask_per`.
11. **GetFinalTaxonomy** — assemble the per-element taxonomy table, joining the
    contributing accessions with ` AND ` when a merged element is supported by
    more than one reference protein.
12. **TagElements** — flag elements that overlap an element of a **different**
    family (within 100 nt) as `overlaped`, others as `unique`, and compute
    `Average_pident`.
13. **GetLength + GetBed + BedFlank + GetFasta** — extract `--flank` nt on each
    side of every element with `bedtools slop`.

## Options reference

### Required inputs

| Option | Meaning |
|--------|---------|
| `-in/--genome_file` | Input genome FASTA (nucleotides). |
| `-od/--outdir` | Output directory (created if missing). |
| `-db/--database` | Reference protein database FASTA. This defines what the run searches for — commonly viral or bacterial proteins, but any donor lineage works. |
| `-mt/--dbmetadata` | Protein metadata CSV for `-db`; supplies every taxonomic name in the output — see [the required format](databases.md#the-metadata-csv-format). |
| `-bt/--hostgenesbaits` | Host-gene bait proteins FASTA — the background candidates are filtered against. |

### Search & filtering

| Option | Default | Meaning |
|--------|---------|---------|
| `-md/--mode` | `blastx` | `blastx` or a DIAMOND sensitivity (`fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive`). DIAMOND is faster but less sensitive — see [The sensitivity trade-off](#diamond-sensitivity). |
| `-ln/--length` | `10000` | Minimum contig length for the search. |
| `-rj/--range_junction` | `100` | Range used to collapse redundant hits — see [Custom arguments](#range-junction-arg). |
| `-p/--threads` | `1` | Threads for the search steps. |
| `-id/--index_databases` | off | Build the BLAST/DIAMOND indexes for `-db` and `-bt`. |

### Element assembly

| Option | Default | Meaning |
|--------|---------|---------|
| `-lm/--limit` | `1` | Bases used to merge neighbouring same-taxon elements (`bedtools merge -d`). |
| `-ml/--merge_level` | `genus` | Taxonomic level (`family`/`genus`) used to decide "same taxon". |
| `-fl/--flank` | `10000` | Flanking-region length extracted on each side. |

### Masking & output

| Option | Default | Meaning |
|--------|---------|---------|
| `-mp/--mask_per` | `50` | Lowercase/ambiguous-base percentage above which a region counts as repetitive. |
| `-cm/--clean_masked` | off | Also emit mask-cleaned outputs (`*.cleaned.*`). |
| `-pr/--prefix` | input filename | Prefix for output files and element names. |
| `-rm/--removetmp` | off | Delete intermediates instead of archiving them under `tmp_files/`. |

See [Outputs](output.md) for the files produced and
[Custom arguments](custom-arguments.md) for the merge and junction behaviour with
worked examples.
