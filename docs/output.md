# Outputs

A `screening` run writes its main results directly into `--outdir`, plus a JSON
run log and (unless `--removetmp`) an archive of intermediates.

## Main outputs

| File | Contents |
|------|----------|
| `PREFIX.EEs.fa` | Endogenous element nucleotide sequences. |
| `PREFIX.EEs.tax.tsv` | Endogenous element taxonomy table (one row per element). |
| `PREFIX.EEs.gff3` | Endogenous element annotation in GFF3 format. |
| `PREFIX.EEs.flanks.fa` | Endogenous elements plus flanking regions (`--flank` nt each side). |
| `eefinder.log` | JSON run summary (see below). |

With `--clean_masked`, the mask-cleaned equivalents are also written:
`PREFIX.EEs.cleaned.fa` and `PREFIX.EEs.cleaned.tax.tsv`.

With `--overlap longest|targets`, the elements filtered out of the final results
are preserved under `tmp_outputs/` as `PREFIX.EEs.removed.{fa,tax.tsv}`.

## The taxonomy table (`PREFIX.EEs.tax.tsv`)

| Column | Meaning |
|--------|---------|
| `Element-ID` | `CONTIG:START-END` identifier of the element (BED-style, 0-based half-open). The run prefix is **not** part of it — see the note below. |
| `Sense` | Strand of the element, `pos` or `neg`. |
| `Protein-IDs` | Accession(s) of the best-matching reference protein(s). |
| `Protein-Products` | Product name(s) from the metadata CSV. |
| `Molecule_type` | Genome composition of the source (e.g. `ssRNA(+)`). |
| `Family` / `Genus` / `Species` | Taxonomic assignment from the metadata CSV. |
| `Host` | Host recorded for the reference. |
| `Overlaped_Element_ID` | Element(s) this one overlaps, if any. |
| `tag` | `overlaped` or `unique`. |
| `Average_pident` | Mean percent identity of the hits backing the element. |

```{note}
The three places an element is named use two spellings. The taxonomy table's
`Element-ID` drops the run prefix (`ctg_1913:1754-2689`), while the `EEs.fa` /
`EEs.flanks.fa` headers and the GFF3 `ID` attribute carry it
(`Ae_aeg_Aag2_ctg_1913/ctg_1913:1754-2689`), so the sequences cross-reference the
genome the run was given. Prepend `{prefix}/` to an `Element-ID` to match a
FASTA record.
```

## The GFF3 annotation (`PREFIX.EEs.gff3`)

Each row of the taxonomy table becomes one feature, sorted by sequence id then
start. The `Element-ID` supplies the sequence id and the coordinates, converted
from the BED-style 0-based half-open range EEfinder carries internally to GFF3's
1-based inclusive convention. `Average_pident` becomes the feature score, and
`--analysis` selects the feature type (`endogenous_viral_element` or
`endogenous_bacterial_element`).

Column 9 carries the taxonomy as attributes — reserved tags capitalised per the
spec, custom tags lower-case:

| Attribute | Source column |
|-----------|---------------|
| `ID` | `Element-ID`, prefixed with the run prefix |
| `Name`, `species` | `Species` |
| `family` / `genus` | `Family` / `Genus` |
| `molecule_type` | `Molecule_type` |
| `product` | `Protein-Products` |
| `protein_ids` | `Protein-IDs` |
| `host` | `Host` |
| `overlap_status` | `tag` (`overlaped` / `unique`) |

GFF3-reserved characters (`;`, `=`, `&`, `,`, tabs, newlines) are
percent-encoded in attribute values.

## The cleaned outputs (`--clean_masked`)

`PREFIX.EEs.cleaned.fa` holds the subset of elements whose lower-case (soft-masked)
plus `N` content is at or below `--mask_per` percent, and
`PREFIX.EEs.cleaned.tax.tsv` is the taxonomy table restricted to those elements.
Both are a **subset** of the unfiltered pair, which is always written — the flag
adds a view rather than replacing one.

## The run log (`eefinder.log`)

`eefinder.log` is a JSON document recording:

- `eefinder_version` — the installed EEfinder version;
- `arguments` — the resolved run arguments (including `translation_method`);
- `dependencies` — detected versions of bedtools, BLAST, DIAMOND, python, numpy
  and pandas, each flagged if it differs from the `env.yml` pin;
- per-step and total timing information.

```{tip}
When EEfinder is installed outside its source tree, point the dependency-drift
check at the reference file with `export EEFINDER_ENV_YML=/path/to/env.yml`.
```

## Intermediate files (`tmp_files/`)

Unless `--removetmp` is given, the intermediates are archived under `tmp_files/`.
Their names accrete suffixes as they pass through the pipeline, so you can trace
exactly which step produced each file. The prefixing and length filtering happen
in one pass, so there is a single `PREFIX.rn.fmt` and no intermediate `PREFIX.rn`:

```text
outdir/
├── eefinder.log
├── PREFIX.EEs.fa
├── PREFIX.EEs.tax.tsv
├── PREFIX.EEs.gff3
├── PREFIX.EEs.flanks.fa
└── tmp_files/
    ├── PREFIX.rn.fmt                 # prefixed headers + length-filtered
    ├── PREFIX.rn.fmt.blastx          # similarity search (main)
    ├── PREFIX.rn.fmt.blastx.filtred  # redundant-hit filter
    ├── PREFIX.rn.fmt.blastx.filtred.bed
    ├── PREFIX.rn.fmt.blastx.filtred.bed.fasta          # putative EEs
    ├── PREFIX.rn.fmt.blastx.filtred.bed.fasta.blastx   # host-bait search
    ├── ...
    └── PREFIX.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge...
```

With the prediction-based translation methods (`gv`/`rv`/`gv-rv`), the
predicted-protein coordinates TSVs also appear here (e.g.
`PREFIX.rn.fmt.pred.coords.tsv` for the main search and
`PREFIX.rn.fmt.blastx.filtred.bed.fasta.pred.coords.tsv` for the host-bait
search) — evidence that the translation method was applied to both searches.
