# Outputs

A run writes its main results directly into `--outdir`, plus a JSON run log and
(unless `--removetmp`) an archive of intermediates.

## Main outputs

| File | Contents |
|------|----------|
| `PREFIX.EEs.fa` | Endogenous element nucleotide sequences. |
| `PREFIX.EEs.tax.tsv` | Endogenous element taxonomy table (one row per element). |
| `PREFIX.EEs.flanks.fa` | Endogenous elements plus flanking regions (`--flank` nt each side). |
| `eefinder.log` | JSON run summary (see below). |

With `--clean_masked`, the mask-cleaned equivalents are written alongside them:
`PREFIX.EEs.cleaned.fa` and `PREFIX.EEs.cleaned.tax.tsv`.

```text
outdir/
├── eefinder.log
├── PREFIX.EEs.fa
├── PREFIX.EEs.tax.tsv
├── PREFIX.EEs.flanks.fa
├── PREFIX.EEs.cleaned.fa          # only with --clean_masked
├── PREFIX.EEs.cleaned.tax.tsv     # only with --clean_masked
└── tmp_files/                     # unless --removetmp
```

## The taxonomy table (`PREFIX.EEs.tax.tsv`)

| Column | Meaning |
|--------|---------|
| `Element-ID` | `CONTIG:START-END` identifier of the element (BED-style, 0-based half-open). The run prefix is **not** part of it — see the note below. |
| `Sense` | Strand of the element, `pos` or `neg`. |
| `Protein-IDs` | Accession(s) of the supporting reference protein(s), each followed by `\|` and the percent identity of its hit. |
| `Protein-Products` | Product name(s) from the metadata CSV. |
| `Molecule_type` | Genome composition of the source (e.g. `ssRNA(+)`, `dsDNA-RT`). |
| `Family` / `Genus` / `Species` | Taxonomic assignment from the metadata CSV. |
| `Host` | Host recorded for the reference. |
| `Overlaped_Element_ID` | Comma-separated element(s) of a **different** family this one overlaps (within 100 nt), if any. |
| `tag` | `overlaped` when `Overlaped_Element_ID` is non-empty, otherwise `unique`. |
| `Average_pident` | Mean of the percent identities listed in `Protein-IDs`, rounded to one decimal. |

An element merged from several fragments carries all its supporting references:
`Protein-IDs` entries are separated by ` | `, while the parallel taxonomy columns
(`Protein-Products`, `Genus`, `Species`, `Host`) join their values with ` AND `.

```text
Element-ID           Sense  Protein-IDs                                Family          Genus                            tag        Average_pident
ctg_1913:1754-2689   neg    YP_006732334.1|30.573 | NP_068729.1|28.239 Caulimoviridae  Soymovirus AND Caulimovirus      overlaped  29.4
ctg_1913:70282-70728 neg    YP_009666257.1|30.719                      Chuviridae      Scarabeuvirus                    unique     30.7
```

Missing metadata is filled rather than left blank: an unresolved `Family`,
`Genus` or `Species` becomes `Unclassified`, and an unresolved `Host` becomes
`Undefined`.

```{note}
The two places an element is named use two spellings. The taxonomy table's
`Element-ID` drops the run prefix (`ctg_1913:1754-2689`), while the `EEs.fa` /
`EEs.flanks.fa` headers carry it
(`Ae_aeg_Aag2_ctg_1913/ctg_1913:1754-2689`), so the sequences cross-reference the
genome the run was given. Prepend `{prefix}/` to an `Element-ID` to match a
FASTA record.
```

## The cleaned outputs (`--clean_masked`)

`PREFIX.EEs.cleaned.fa` holds the subset of elements whose lower-case (soft-masked)
plus `N` content is at or below `--mask_per` percent, and
`PREFIX.EEs.cleaned.tax.tsv` is the taxonomy table restricted to those elements.
Both are a **subset** of the unfiltered pair, which is always written — the flag
adds a view rather than replacing one.

## The run log (`eefinder.log`)

`eefinder.log` is a JSON document recording:

- `arguments` — the resolved run arguments;
- `start_time`, `end_time`, `total_time_minutes` — wall-clock timing of the run;
- `steps_information` — one entry per pipeline step, each with its own start/end
  time, duration and a human-readable message describing what it did with which
  parameters.

```json
{
    "arguments": {
        "genome_file": "test_files/Ae_aeg_Aag2_ctg_1913.fasta",
        "prefix": "Ae_aeg_Aag2_ctg_1913",
        "outdir": "results_test",
        "mode": "blastx",
        "length": 1000,
        "...": "..."
    },
    "start_time": "2024-10-01 12:00:00",
    "total_time_minutes": "1.8421",
    "steps_information": [
        {
            "step": "Prepare input data",
            "total_time_minutes": "0.0032",
            "message": "Ae_aeg_Aag2_ctg_1913 prefix included in ... sequences header ..."
        }
    ]
}
```

## Intermediate files (`tmp_files/`)

Unless `--removetmp` is given, the intermediates are archived under `tmp_files/`.
Their names accrete suffixes as they pass through the pipeline, so you can trace
exactly which step produced each file:

```text
tmp_files/
├── PREFIX.rn                        # headers prefixed with PREFIX/
├── PREFIX.rn.fmt                    # contigs below --length removed
├── PREFIX.rn.fmt.fai                # samtools index written by bedtools getfasta
├── PREFIX.rn.fmt.rn.fmt.lenght      # contig lengths, for bedtools slop
├── PREFIX.rn.fmt.blastx             # main similarity search
├── PREFIX.rn.fmt.blastx.csv         # the same hits with the working columns added
├── PREFIX.rn.fmt.blastx.filtred     # redundant-hit filter (--range_junction)
├── PREFIX.rn.fmt.blastx.filtred.bed
├── PREFIX.rn.fmt.blastx.filtred.bed.fasta                    # putative EEs
├── PREFIX.rn.fmt.blastx.filtred.bed.fasta.blastx             # host-bait search
├── PREFIX.rn.fmt.blastx.filtred.bed.fasta.blastx.csv
├── PREFIX.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred
├── ....filtred.concat               # EE + host hits, ranked by bitscore
├── ....concat.nr                    # best hit per element, host winners dropped
├── ....concat.nr.tax                # joined with the -mt metadata
├── ....tax.bed                      # annotated with the --merge_level taxon
├── ....tax.bed.merge                # bedtools merge -d --limit
├── ....merge.fmt                    # annotation stripped back to coordinates
├── ....merge.fmt.fa.bed
└── ....merge.fmt.fa.bed.flank       # coordinates grown by --flank
```

The three merged-element files that become the main outputs
(`....merge.fmt.fa`, `....merge.fmt.fa.tax` and `....merge.fmt.fa.bed.flank.fasta`)
are renamed into `PREFIX.EEs.fa`, `PREFIX.EEs.tax.tsv` and
`PREFIX.EEs.flanks.fa` rather than archived.
