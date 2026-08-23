# Changelog

All notable changes to EEfinder are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/).


## [2.0.0]

Major release. The CLI became a command group, a database-acquisition command
was added, the similarity search gained alternative translation methods and
overlap resolution, and the project gained GFF3 output, a run/dependency audit,
a full test suite, a Read-the-Docs docs and a `pyproject.toml` build.

### Added
- **`get-databases` command** (`eefinder/get_databases.py`) — downloads the
  RefSeq protein databases EEfinder needs via the NCBI `datasets` CLI
  (`ncbi-datasets-cli>=18.1`, pinned in `env.yml`; earlier builds fail on viral
  downloads with `ACELLULAR_ROOT is not a valid V2reportsRankType`). It is a
  command group with one subcommand per database:
  - `virus` (taxon default `10239`) and `bacteria` (taxon default `2`) write a
    protein FASTA **and** a metadata CSV
    (`Accession,Species,Genus,Family,Molecule_type,Protein,Host`);
  - `host` (taxon required) writes the `-bt` baits FASTA only.

  The CSV is rebuilt from the downloaded `protein.faa` headers joined with the
  `data_report.jsonl` taxonomy report: `Host` from the report, `Genus`/`Family`
  inferred from the report's unranked lineage by ICTV name suffix
  (`-viridae` / single-word `-virus`), and `Molecule_type` — absent from the
  report — filled from a bundled
  [ICTV genome-composition table](https://ictv.global/virus-properties)
  (`eefinder/data/ictv_genome_composition.tsv`). Options: `--exclude-uninformative`
  (drop `hypothetical`/`uncharacterized protein` records), `--standardize-proteins`
  (below), `--cluster` (below) and `--refseq`. Each run writes a JSON log
  `{outdir}/{prefix}.log` with a `sequence_counts` block
  (`downloaded`/`excluded_uninformative`/`clustered_identical`/
  `dropped_standardization`/`kept`). Pure helpers are unit-tested with no
  network / `datasets` binary.
- **`get-databases --cluster/--no-cluster` (on by default).** Collapses
  100%-identical / 100%-coverage duplicate proteins with `cd-hit`
  (`-c 1.0 -aL 1.0 -aS 1.0`) before building the metadata CSV — a lossless dedup
  that shrinks the database and speeds up the search. Needs `cd-hit` on `PATH`
  (checked up front, with a clear error pointing to `--no-cluster`).
- **`get-databases --standardize-proteins` protein-name standardization**
  (`eefinder/normalization.py`, `standardize_protein(name, mol_type, target)`).
  Dispatches per target: `virus` applies a bundled first-draft canonical-name map
  (`eefinder/data/viral_proteins.tsv`), collapsing synonyms per `Molecule_type`
  scope (e.g. every RdRp spelling — including compound names like `P2-RdRp` or
  `CP/RdRp fusion` — → `RdRp`; the Capsid, Glycoprotein and Nucleocapsid variants
  to their canonical forms); `bacteria`/`host` get generic cleaning only. All
  targets share a cleaning pass: leaked NCBI `[key=value]` tags
  (`[organism=…]`, `[gbkey=CDS]`, one level of nesting) and the characters
  `:,/\?!` + quotes are removed, molecular-weight tokens are normalised
  (`33 kDa`, `33-kDa`, `33K-like protein` → `33 kDa protein`), misspellings are
  fixed so variants converge, leading hedging qualifiers and `CDS:`/`ORF:`
  directives are stripped, and the first letter is capitalised. Records that
  reduce to `Unknown` (bare `CDS`/`ORF`) or begin with `hypothetical` are dropped
  from both the CSV and the FASTA (kept in sync via `filter_fasta_by_ids`).
- **`screening --translation_method {default,gv,rv,gv-rv}` (`-tm`).** Controls
  how proteins are obtained for the similarity search, applied consistently to
  **both** searches in a run (the main EE search and the host-bait search) via a
  single value threaded into `SimilaritySearch`, so they can never diverge.
  `default` keeps the six-frame `blastx`/`diamond blastx`; `gv`/`rv` predict
  proteins from the nucleotide query with
  [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) /
  [pyrodigal-rv](https://github.com/LanderDC/pyrodigal-rv) then align with
  `blastp`/`diamond blastp`; `gv-rv` runs both predictors and drops redundancy
  with `cd-hit` (100%/100%). New `eefinder/translation.py` handles prediction
  (writing a coordinates TSV `protein_id,contig,start,end,strand,tool`) and a
  **coordinate traceback** that maps each protein hit's amino-acid span back to
  nucleotide coordinates on the source contig — so `SimilaritySearch` always
  emits the same `{query}.blastx` schema and every downstream step is unchanged.
  New deps: `pyrodigal-gv`, `pyrodigal-rv` (pip) and `cd-hit` (conda).
- **`screening --overlap {keep,longest,targets}` (`-ov`)** — resolution strategy
  for elements tagged `overlap_status=overlaped` (`eefinder/overlap.py`):
  `keep` (default) keeps everything; `longest` keeps the longest of each overlap
  cluster; `targets` resolves each cluster by **exactly one** of two repeatable,
  mutually-exclusive family lists — `--target_families`/`-tf` (keep-list) or
  `--non_target_families`/`-ntf` (drop-list, never wipes an all-listed cluster).
  Filtering applies to every final result; the filtered-out elements are
  preserved under `tmp_outputs/` (`{prefix}.EEs.removed.*`).
- **GFF3 output** (`eefinder/gff.py`, `WriteGFF3`). The pipeline emits
  `{prefix}.EEs.gff3` (and `{prefix}.EEs.cleaned.gff3` under `--clean_masked`),
  converting the taxonomy table to GFF3 (BED 0-based → 1-based, `Sense` → strand,
  `Average_pident` → score, taxonomy fields as escaped column-9 attributes, `ID`
  = `{prefix}/{Element-ID}` to match the FASTA headers).
- **`screening --analysis {virus,bacteria}` (`-an`)** selects the GFF3 feature
  type (`endogenous_viral_element` / `endogenous_bacterial_element`).
- **Run-context / dependency audit in `eefinder.log`** (`eefinder/versions.py`).
  Detects the versions of the external tools (bedtools, BLAST, DIAMOND) and
  Python libraries used and compares them to the `env.yml` pins (found via the
  package location or `EEFINDER_ENV_YML`). The log records `eefinder_version`, a
  `dependencies` list (`ok`/`mismatch`/`not-found`/`unpinned`) and the host
  `system` context; a startup header warns about any drift.
- **`--debug` flag** on `screening` and every `get-databases` subcommand — lowers
  the `eefinder` logger to DEBUG (`log.enable_debug()`) for verbose traces. Off
  by default.
- **`pytest` test suite** under `tests/` — unit tests for every data-processing
  step (synthetic inputs in `tmp_path`, no binaries) plus a scenario-driven set
  of end-to-end integration tests against `test_files/` (auto-skipped without the
  external binaries), including byte-for-byte golden comparisons of the main
  outputs under `test_files/expected_results/{default,gv-rv}/`. `pytest
  --update-test` regenerates the golden files.
- **GitHub Actions workflow** (`.github/workflows/tests.yml`) running the
  binary-free unit tests (`pytest -m "not integration"`) on every pull request to
  `master`, installing only the pip deps (no BLAST/DIAMOND/bedtools).
- **Read-the-Docs documentation site** under `docs/` (Sphinx + MyST +
  `sphinx_rtd_theme`, `.readthedocs.yaml`) with pages for installation,
  `get-databases`, running the pipeline, translation methods, overlap resolution,
  custom arguments, outputs, testing and a developer guide (builds warning-clean).
- **`CLAUDE.md`** project guidance and **`CHANGELOG.md`** (this file).

### Changed
- **The CLI is now a command group; the pipeline moved under `screening`.** The
  console entry point is the `cli` group, so the pipeline is invoked as
  **`eefinder screening <options>`** instead of `eefinder <options>`. Options are
  otherwise unchanged. **(Breaking.)**
- **`--merge_level` default changed from `genus` to `family`**; the default GFF3
  feature type is `endogenous_viral_element` (was `translated_nucleotide_match`).
- **Build system migrated to `pyproject.toml`** (hatchling backend, `[project]`
  table with runtime deps + `dev` extra + the `eefinder` console script);
  `setup.py`/`MANIFEST.in` removed. `pyproject.toml` also holds the black + pytest
  config. `pip install .` now pulls the Python runtime deps (click, biopython,
  pandas<2, numpy<2); the external binaries still come from `env.yml`.
- **Updated pinned tool versions in `env.yml`**: bedtools 2.27.1 → 2.31.1, BLAST
  2.5.0 → 2.17.0, DIAMOND 2.0.15 → 2.2.3, plus `ncbi-datasets-cli`, `cd-hit`,
  `pyrodigal-gv` and `pyrodigal-rv`. `env.yml` now pins only the direct
  dependencies and lets conda resolve the rest (the old frozen transitive pins
  conflicted with BLAST 2.17's openssl 3.x). BLAST 2.17.0 is unavailable on
  osx-arm64 — pin `blast=2.16.0` there.
- **`eefinder/__init__.py` resolves `__version__` via `importlib.metadata`**
  instead of `pkg_resources`, removing the deprecation warning and the implicit
  `setuptools` runtime pin.
- **Refactored the pipeline for clean code, dataclasses, NumPy-style docstrings
  and type hints throughout**, without changing default-path outputs (verified
  byte-identical on `test_files/`): the run-log dicts became the `StepInfo` /
  `RunArguments` / `RunInfo` dataclasses; shared constants were extracted; dead
  code removed; kept parseable under the pinned Python 3.9.
- `.gitignore` covers the full set of BLAST/DIAMOND index suffixes plus Python
  build/pytest artifacts. Subprocess command strings reformatted one parameter
  per line.

### Removed
- **`setup.py` and `MANIFEST.in`** — build metadata moved to `pyproject.toml`.
- **The top-level pipeline invocation `eefinder <options>`** — use
  `eefinder screening <options>`. **(Breaking.)**

### Fixed
- **`--clean_masked` produced an empty cleaned taxonomy table.** Cleaned FASTA
  record IDs keep the `PREFIX/` that `TagElements` strips from the taxonomy
  `Element-ID`, so the id comparison never matched. IDs are now compared with the
  prefix removed (and the source table's 12-column header preserved), so
  `*.EEs.cleaned.tax.tsv` is populated.
- **`RunInfo.merge_level` was mis-mapped** to the `length` argument in the run
  log; it now reports the real `--merge_level`.
- Corrected the run-log timestamp format string (stray `%` in `%H:%M%:%S`).

### Notes
- The default `blastx` mode is the tested, reliable path. The DIAMOND modes can
  fail silently on newer `diamond` builds (subprocess stderr is discarded); verify
  the `diamond` build pinned in `env.yml` if a DIAMOND run produces no hits.

[2.0.0]: https://github.com/WallauBioinfo/EEfinder/compare/v1.1.1...v2.0.0
[1.1.1]: https://github.com/WallauBioinfo/EEfinder/releases/tag/v1.1.1
