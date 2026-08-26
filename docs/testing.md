# Testing EEfinder

EEfinder ships with a [pytest](https://docs.pytest.org/) suite under `tests/`.
It has two layers:

| Layer | Marker | Needs BLAST/DIAMOND/bedtools? | What it covers |
|-------|--------|-------------------------------|----------------|
| Unit  | (none) | No | The pure data-processing steps, with small synthetic inputs written to `tmp_path`. |
| Integration | `integration` | Yes | End-to-end runs of the `eefinder` CLI against the example files in `test_files/`. |

## Install the development dependencies

The unit tests need only the runtime Python dependencies plus pytest. The
integration tests additionally need the external binaries, which come from the
conda environment.

```bash
# runtime env + external binaries (blast, diamond, bedtools)
micromamba env create -f env.yml        # or: conda env create -f env.yml
micromamba activate EEfinder
pip install -e .

# test + lint tooling
pip install -r requirements-dev.txt
```

`requirements-dev.txt` pins:

- `pytest` — test runner
- `black` — code formatter / linter

## Running the tests

```bash
# everything
pytest

# unit tests only — no external binaries required
pytest -m "not integration"

# integration (end-to-end) tests only
pytest -m integration
```

The integration tests call `shutil.which(...)` for `eefinder`, `blastx`,
`makeblastdb` and `bedtools`; if any is missing they are **skipped** (not
failed), so `pytest -m "not integration"` runs cleanly on a bare Python install.

## What the tests check

### Unit tests

- **`test_prepare_data.py`** — `InsertPrefix` rewrites every FASTA header to
  `>PREFIX/…` and leaves sequences untouched; `PrepareGenome` does the prefixing
  and the length cutoff in one pass, byte-for-byte matching the two-step result.
- **`test_clean_data.py`** — `RemoveShortSequences` applies the length cutoff
  (inclusive); `MaskClean` drops soft-masked sequences at the given threshold.
- **`test_filter_table.py`** — `FilterTable` collapses redundant hits, swaps
  coordinates for negative-strand hits, drops hits shorter than 33 aa, and emits
  the `.filtred`/`.filtred.bed` files with the expected `bed_name`/`tag`.
- **`test_taxonomy.py`** — (the taxonomy *table*, not `taxon_exclusion.py`)
  `GetTaxonomy` left-joins the metadata CSV;
  `CompareResults` removes putative EEs that hit the host baits with a higher
  bitscore.
- **`test_bed.py`** — `GetBed`, `RemoveAnnotation` and `GetAnnotBed` produce the
  expected BED coordinates and annotation strings.
- **`test_tag_elements.py`** — `TagElements` flags overlapping vs unique
  elements and computes `Average_pident`.
- **`test_gff.py`** — `WriteGFF3` emits a valid GFF3 (1-based coordinates,
  strand mapping, percent-escaped attributes, custom source/type).
- **`test_versions.py`** — parsing `env.yml` pins and classifying detected vs
  expected dependency versions (`ok`/`mismatch`/`not-found`/`unpinned`).
- **`test_get_length.py`** — `GetLength` writes `id<TAB>length` per record.
- **`test_get_databases.py`** — the `get-databases` helpers: protein-header and
  `data_report.jsonl` parsing, metadata-CSV assembly, `datasets` command
  construction, `protein.faa` concatenation, and the CLI's missing-`datasets`
  guard (no network or `datasets` binary required).
- **`test_overlap.py`** — `elements_to_remove` for the `keep`/`longest`/
  `targets` strategies (per-cluster keep-list and drop-list resolution) and
  `FilterOverlap` splitting kept vs removed records.
- **`test_taxon_exclusion.py`** — taxonomy expansion for `--exclude-taxon`:
  sibling collection along the lineage, composing several exclusions, exclusions
  outside the requested taxon, `--inputfile` batching, and the error cases. The
  `datasets` CLI is faked, so no network is used.
- **`test_translation.py`** — protein prediction, `cd-hit` clustering and the
  amino-acid → nucleotide coordinate traceback.
- **`test_progress.py`** / **`test_download_hang.py`** — the progress display,
  retry/backoff and stall detection, including the cases where `datasets`
  delivers a complete package but never exits.
- **`test_utils.py`** — `check_outdir`, `step_info`, `running_info` helpers.

### Integration tests

`test_integration.py` runs the CLI on `test_files/` (base command: `blastx`,
`-ln 1000 -id -p 2 -lm 100`). It is a small, scenario-driven set rather than an
exhaustive parameter sweep — each test covers one behaviour worth protecting:

- **`test_matches_expected_results`** — the documented `default` run reproduces
  the four main outputs (`*.EEs.fa`, `*.EEs.tax.tsv`, `*.EEs.flanks.fa`,
  `*.EEs.gff3`) **byte-for-byte** against the golden copies in
  [`test_files/expected_results/default/`](https://github.com/WallauBioinfo/EEfinder/tree/master/test_files/expected_results/default)
  (the timestamped `eefinder.log` and temporary files are ignored).
- **`test_matches_expected_results_gv_rv`** — the same byte-for-byte comparison
  for a full `--translation_method gv-rv` run against
  [`test_files/expected_results/gv-rv/`](https://github.com/WallauBioinfo/EEfinder/tree/master/test_files/expected_results/gv-rv);
  gated on `pyrodigal-gv`/`-rv` + `blastp` + `cd-hit` and skipped otherwise.
- **`test_translation_method_gv_drives_both_searches`** — `-tm gv` predicts
  proteins for **both** the main and the host-bait search (a coordinates TSV
  exists for each), with the schema unchanged.
- **`test_clean_masked_is_subset_of_full_run`** — `--clean_masked` produces a
  populated cleaned table that is a subset of the full run.
- **`test_merge_limit_controls_merging`** — a larger `--limit` merges
  neighbouring same-taxon elements, reducing the element count.
- **`test_family_merge_level_runs_end_to_end`** — the `--merge_level family`
  branch runs and yields a valid taxonomy table.
- **`test_overlap_longest_filters_and_preserves_removed`** — `--overlap longest`
  drops overlaping elements from the final results while preserving them under
  `tmp_outputs/`; kept and removed elements partition the unfiltered run.
- **`test_overlap_targets_requires_exactly_one_family_list`** — `--overlap
  targets` with neither (or both) of `--target_families` / `--non_target_families`
  exits non-zero with a helpful message.
- **`test_overlap_targets_non_target_families_drops_that_family`** — `--overlap
  targets --non_target_families` removes the listed family from the results while
  preserving it under `tmp_outputs/`.
- **`test_removetmp_removes_intermediate_files`** — `--removetmp` deletes the
  intermediates instead of archiving them.

If you deliberately change the pipeline output (for example after bumping a
dependency version), refresh the golden files instead of editing them by hand:

```bash
pytest -m integration --update-test
```

`--update-test` makes the golden tests overwrite their respective
`test_files/expected_results/{default,gv-rv}/` directory with the freshly
produced outputs and skip the comparison. Review the diff and commit the
refreshed files.

## Formatting

Format the code (and check in CI) with black:

```bash
black eefinder tests      # apply
black --check eefinder tests   # verify only
```
