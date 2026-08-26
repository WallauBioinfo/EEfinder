# Running the pipeline (`screening`)

`screening` is the EE-finding pipeline. It takes a genome FASTA plus the protein
database, its metadata CSV, and the host-gene baits, and produces the endogenous
element sequences, taxonomy table, GFF3 annotation, and flanking regions.

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
eefinder screening \
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
- `-id` builds the BLAST/DIAMOND indexes for the databases (needed on the first
  run against a given database).
- `-p 2` uses two threads.
- `-lm 100` merges same-taxon elements within 100 nt.

## Pipeline steps

The `screening` command orchestrates these steps (each a small side-effect class;
files flow through disk with accreting suffixes `.rn.fmt`, `.blastx`,
`.filtred`, `.bed`, `.tax`, …):

1. **PrepareGenome** — prefix every FASTA header (`>PREFIX/…`) **and** drop
   contigs below `--length`, in a single pass writing only `{prefix}.rn.fmt`.
   (`InsertPrefix` and `RemoveShortSequences` still exist as standalone classes;
   chaining them wrote the genome to disk twice.)
2. **MakeDB** — build BLAST or DIAMOND databases (`--index_databases`).
3. **SimilaritySearch** — the similarity search, run **twice** (main EE search +
   host-bait search). `--translation_method` controls both — see
   [Translation methods](translation-methods.md).
4. **FilterTable** — filter redundant hits by `qseqid`/range/sense.
5. **GetFasta** — extract putative EE sequences (bedtools).
6. **CompareResults** — drop EEs that hit host baits harder.
7. **GetTaxonomy** — join hits to the metadata CSV, build the taxonomy table.
8. **MergeBed** — merge truncated elements of the same genus/family
   (`--merge_level`).
9. **MaskClean** — optional soft-mask filter (`--clean_masked`).
10. **TagElements** — flag overlapping elements, add `Average_pident`.
11. **FilterOverlap** — resolve overlaps by the chosen strategy — see
    [Overlap resolution](overlap.md).
12. **WriteGFF3** — write the EE taxonomy table as a GFF3 annotation.
13. **GetLength + BedFlank + GetFasta** — extract flanking regions (`--flank`).

(diamond-sensitivity)=
## The sensitivity trade-off

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

## Options reference

### Required inputs

| Option | Meaning |
|--------|---------|
| `-in/--genome_file` | Input genome FASTA (nucleotides). |
| `-od/--outdir` | Output directory. |
| `-db/--database` | Protein database FASTA (virus or bacteria). |
| `-mt/--dbmetadata` | Protein metadata CSV for `-db`. |
| `-bt/--hostgenesbaits` | Host-gene bait proteins FASTA. |

### Search & filtering

| Option | Default | Meaning |
|--------|---------|---------|
| `-md/--mode` | `blastx` | `blastx` or a DIAMOND sensitivity (`fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive`). DIAMOND is faster but less sensitive — see [The sensitivity trade-off](#diamond-sensitivity). |
| `-tm/--translation_method` | `default` | `default`/`gv`/`rv`/`gv-rv` — see [Translation methods](translation-methods.md). |
| `-ln/--length` | `10000` | Minimum contig length for the search. |
| `-rj/--range_junction` | `100` | Range for junction of redundant hits — see [Custom arguments](custom-arguments.md). |
| `-p/--threads` | `1` | Threads for multi-threaded steps. |
| `-id/--index_databases` | off | Build the BLAST/DIAMOND indexes for the databases. |

### Element assembly

| Option | Default | Meaning |
|--------|---------|---------|
| `-lm/--limit` | `1` | Bases used to merge neighbouring same-taxon elements (bedtools merge). |
| `-ml/--merge_level` | `family` | Taxonomic level (`family`/`genus`) to merge by. |
| `-fl/--flank` | `10000` | Flanking-region length to extract. |
| `-ov/--overlap` | `keep` | `keep`/`longest`/`targets` — see [Overlap resolution](overlap.md). |
| `-tf/--target_families` | — | Family to KEEP (repeatable) with `--overlap targets`. |
| `-ntf/--non_target_families` | — | Family to DROP (repeatable) with `--overlap targets`. |

### Masking & output

| Option | Default | Meaning |
|--------|---------|---------|
| `-mp/--mask_per` | `50` | Lowercase-percentage threshold to call a region repetitive. |
| `-cm/--clean_masked` | off | Also emit mask-cleaned outputs (`*.cleaned.*`). |
| `-an/--analysis` | `virus` | GFF3 feature type (`virus` → `endogenous_viral_element`, `bacteria` → `endogenous_bacterial_element`). |
| `-pr/--prefix` | input filename | Prefix for output files and Element-IDs. |
| `-rm/--removetmp` | off | Delete intermediates instead of archiving them under `tmp_files/`. |
| `--debug` | off | Emit verbose debug logging (intermediate paths, per-step details). |

```{note}
The default `blastx` mode is the tested, reliable path. The DIAMOND modes
depend on the `diamond` build pinned in `env.yml` and can fail silently if that
build misbehaves — verify it if a DIAMOND run produces no hits.
```

See [Outputs](output.md) for the files produced and
[Custom arguments](custom-arguments.md) for the merge/junction behaviour with
worked examples.
