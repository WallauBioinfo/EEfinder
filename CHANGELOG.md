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
- **`--exclude-taxon` on `get-databases`** (`eefinder/taxon_exclusion.py`, new) — leaves
  a branch of the NCBI taxonomy out of a download **without ever requesting it**.
  SARS-CoV-2 is 9,210,670 of the 15,027,633 viral records in GenBank (61%, and
  99.2% of *Coronaviridae*), so `get-databases virus --all-sequences` otherwise
  spends most of its bandwidth and clustering time on one virus. The `datasets`
  CLI has no exclusion flag, so EEfinder inverts the request: it walks the
  lineage from the requested taxon down to the excluded one and keeps every child
  except the one on the path (`expand_taxon_excluding`). For all viruses minus
  SARS-CoV-2 that is 63 taxa, which fit in one `datasets --inputfile` call;
  longer lists are split across downloads and extracted one package per
  subdirectory (`find_data_reports`, `batch_taxa`). The option is repeatable,
  accepts tax ids or scientific names, defaults to `2697049` (SARS-CoV-2) for
  `virus` and to nothing for `bacteria`/`host`, and is switched off with
  `--exclude-taxon none`. Recorded in the run log as `arguments.exclude_taxa`.
  Verified end to end against *Coronaviridae*: the default build has 1,042
  records and no SARS-CoV-2, `--exclude-taxon none` has 1,080, and the two differ
  by exactly the 38 SARS-CoV-2 accessions — SARS-CoV-**1** (Tor2) is retained.
  Known limit: records attached directly to a rank on the path (`"Betacoronavirus
  sp."`) are unreachable this way — 273 of 5.8 million (0.005%).
- **`get-databases` command** (`eefinder/get_databases.py`) — downloads the
  RefSeq protein databases EEfinder needs via the NCBI `datasets` CLI
  (`ncbi-datasets-cli`, pinned to 18.36.0 in `env.yml`; earlier builds fail on viral
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
- **Documentation pages for the 2.0.0 features.** The Read-the-Docs site
  (Sphinx + MyST + `sphinx_rtd_theme`, `.readthedocs.yaml`) arrived in 1.1.2;
  2.0.0 adds `get-databases`, translation methods, overlap resolution and
  testing pages, rewrites the pipeline page around the `screening` command, and
  audits the existing pages against the code. `databases.md` and `running.md`
  are folded into `get-databases.md` and `screening.md`. Builds warning-clean.
- **`CLAUDE.md`** project guidance and **`CHANGELOG.md`** (this file).
- **Progress reporting for `get-databases`** (`eefinder/progress.py`). The
  `datasets` CLI draws its own transfer progress and EEfinder was swallowing it
  with `capture_output=True`; its output is now mirrored to the terminal as it
  arrives, while still being captured so a failure can quote the tool's own
  message. The steps EEfinder performs itself — archive extraction, `protein.faa`
  merging and metadata building — get `click` progress bars. Progress is shown
  only on a terminal (piped runs stay silent, as before) and
  `EEFINDER_NO_PROGRESS=1` disables it.
- **Retries and stall detection for downloads** (`--attempts`, default 3;
  `--stall-timeout`, default 180s). NCBI transfers fail both by reporting an
  error (`Internal error (invalid zip archive). Please try again`) and by hanging
  with the connection open, which no exit status ever reports. A transfer counts
  as stalled only when **neither** new output from `datasets` **nor** growth of
  the archive is seen within the timeout, so a slow-but-working transfer is not
  mistaken for a hung one; the attempt is then killed, the partial archive
  discarded and the download retried with a growing backoff. Permanent failures
  (a misspelled taxon, an unknown flag) are detected and **not** retried. Both
  settings are recorded in the `{prefix}.log` arguments block.

- **`get-databases --released-before YYYY-MM-DD`** restricts a build to data
  released on or before a date, so a database can be rebuilt as it stood then;
  the date is recorded in `{prefix}.log`. For `bacteria`/`host` the cutoff is the
  NCBI client's own flag (exact, per assembly). Its **virus** subcommand has no
  such flag, so EEfinder applies it from the release dates in the download's
  report — per **organism**, because the protein FASTA carries no link back to
  the record a protein came from (mapping by `proteinCount` was tested and does
  not hold: in *Peribunyaviridae* the totals match but the order does not, in
  *Orthomyxoviridae* the totals disagree). An organism is kept when its earliest
  record predates the cutoff, so a segmented virus whose later segment postdates
  it is not erased. Excluded accessions are recorded in the tracking table with
  `Reason = released_after_cutoff`.
- **`get-databases` writes a `{prefix}.tracking.tsv`** recording the fate of
  every downloaded accession: its product name as delivered and after
  standardisation (with a `Name_changed` flag), whether it was removed and why
  (`uninformative_product`, `identical_duplicate`,
  `product_standardized_to_unknown`), and for clustered sequences the `cd-hit`
  cluster and the representative accession that replaced them. The sequences
  dropped before the FASTA is written appear nowhere else, so this is the only
  record that the download and the database can be reconciled against. The table
  is verified against the finished FASTA, so a record that vanished without any
  step reporting it is still accounted for (`absent_from_final_database`) —
  `cd-hit` silently discards sequences of up to 10 residues (`-l`), which would
  otherwise be reported as kept while missing from the database. Every row also
  carries `Organism_release_date` (the earliest release date among the organism's
  records), so a dated build can be checked rather than trusted.
- **The raw download is deleted once it has been read** — the zip and the
  extracted `*_ncbi/` directory are several GB of duplication for a whole-RefSeq
  run. `--keep-download` keeps them. What a run leaves behind is the FASTA, the
  metadata CSV, the tracking table, the `cd-hit` cluster file
  (`{prefix}.clstr`, new — the members and representative of every cluster) and
  the JSON log.

### Changed
- **The genome is prepared in a single pass** (`prepare_data.PrepareGenome`).
  Prefixing the headers and dropping the short contigs used to be two steps
  chained through a `{prefix}.rn` file, so the whole genome was written to disk
  twice — 3.2 GB of intermediates for a 1.6 GB genome, of which the first copy
  was never read again. Only `{prefix}.rn.fmt` is written now, and the output is
  byte-for-byte what the two steps produced (asserted by a test, and by the
  golden-file integration runs). `InsertPrefix` and `RemoveShortSequences` remain
  for callers that need one operation on its own.
- **The CLI is now a command group; the pipeline moved under `screening`.** The
  console entry point is the `cli` group, so the pipeline is invoked as
  **`eefinder screening <options>`** instead of `eefinder <options>`. Options are
  otherwise unchanged. **(Breaking.)**
- **`--merge_level` default changed from `genus` to `family`.**
- **Packaging.** The console script is repointed at the `cli` group (above);
  1.1.2's packaging metadata (MIT `license`/`license-files`, maintainers, URLs,
  classifiers, `requires-python = ">=3.9,<3.12"`, the sdist manifest) and its
  `publish.yml` Trusted Publishing workflow are otherwise unchanged. The release
  path is verified end to end: `python -m build`, `twine check --strict`,
  install into a clean venv, `eefinder --version` → `eefinder, version 2.0.0`.
- **Dependency bounds corrected.** `biopython` is capped at `<1.86`, which
  **removed** `Bio.Blast.Applications` — still used by `make_database.py` and
  `similarity_analysis.py`. `pandas`/`numpy` were relaxed from `<2` to `<3`
  after confirming the unit suite passes in full against pandas 2.3 / numpy 2.2.
- **`pyproject.toml` extended.** The hatchling build and the `[project]` table
  arrived in 1.1.2; 2.0.0 adds the `dev` extra (pytest + black), the black and
  pytest configuration, and points the console script at the `cli` group instead
  of the `main` command. `pip install .` pulls the Python runtime deps; the
  external binaries still come from `env.yml`.
- **Updated pinned tool versions in `env.yml`**: bedtools 2.27.1 → 2.31.1, BLAST
  2.5.0 → 2.17.0, DIAMOND 2.0.15 → 2.2.3, plus `ncbi-datasets-cli`, `cd-hit`,
  `pyrodigal-gv` and `pyrodigal-rv`. `env.yml` now pins only the direct
  dependencies and lets conda resolve the rest (the old frozen transitive pins
  conflicted with BLAST 2.17's openssl 3.x).
- **Refactored the pipeline for clean code, dataclasses, NumPy-style docstrings
  and type hints throughout**, without changing default-path outputs (verified
  byte-identical on `test_files/`): the run-log dicts became the `StepInfo` /
  `RunArguments` / `RunInfo` dataclasses; shared constants were extracted; dead
  code removed; kept parseable under the pinned Python 3.9.
- `.gitignore` covers the full set of BLAST/DIAMOND index suffixes plus Python
  build/pytest artifacts. Subprocess command strings reformatted one parameter
  per line.

### Removed
- **The top-level pipeline invocation `eefinder <options>`** — the console entry
  point is now the `cli` group, so the pipeline is `eefinder screening
  <options>`. In 1.1.2 the entry point was the single `main` command, which ran
  the pipeline directly. **(Breaking.)**

### Fixed
- **`--clean_masked` produced an empty cleaned taxonomy table.** Cleaned FASTA
  record IDs keep the `PREFIX/` that `TagElements` strips from the taxonomy
  `Element-ID`, so the id comparison never matched. IDs are now compared with the
  prefix removed (and the source table's 12-column header preserved), so
  `*.EEs.cleaned.tax.tsv` is populated.
- **`RunInfo.merge_level` was mis-mapped** to the `length` argument in the run
  log; it now reports the real `--merge_level`.
- Corrected the run-log timestamp format string (stray `%` in `%H:%M%:%S`).
- **`get-databases` hung after the download had finished** — the run sat forever
  on a completed `Validating package files … 100%`. Root cause: the NCBI
  `datasets` client never exits when it inherits an **open standard input**,
  which is what happens in an interactive terminal. Measured on the same
  download, same machine: **3 seconds** with stdin closed, **indefinite** with it
  open. Every external command is now run with stdin closed; nothing here feeds a
  command through stdin, so this costs nothing.
- **The progress display could freeze mid-word** (`Valida…`): the child's pipe
  was read as text in fixed 256-character blocks, which blocks until that many
  characters exist, so a tool redrawing a progress bar appeared stuck until
  enough further output accumulated. The pipe is now read unbuffered and in
  binary, forwarding output the moment it arrives, with incremental decoding so a
  chunk boundary inside a multi-byte character does not corrupt it.
- **Waiting for end-of-file could outlive the command.** A background process a
  command leaves behind holds the write end of the pipe open, so the pipe never
  closes even after the tool has exited. Exit of the command now ends the loop;
  output still in flight is drained first.
- **An HTTP/2 stream reset is retried over HTTP/1.1.** NCBI aborts large
  transfers with `stream error: stream ID 3; INTERNAL_ERROR; received from peer`;
  the same download completes when the Go client is told to avoid HTTP/2
  (`GODEBUG=http2client=0`), so that is what the retry does instead of repeating
  the attempt unchanged. Measured on the full RefSeq viral set: HTTP/1.1 finished
  in 176s where HTTP/2 failed three times in a row.
- **A download that leaves a broken archive is no longer accepted.** NCBI
  occasionally delivers a full-size package that is not a readable zip
  (`118MB invalid zip archive`), sometimes while still reporting success. The
  archive is now verified after the download and a bad one is retried like any
  other failure. Any archive already at the destination is removed **before** the
  first attempt too, so an interrupted run never leaves something to download
  into.
- **A quiet command no longer looks like a freeze.** While a command produces no
  output and its output file does not grow, EEfinder says so every 20 seconds,
  including how long it has been quiet and when it will give up and retry —
  silence was indistinguishable from a hang, so users interrupted runs that were
  about to recover on their own.
- `kill -USR1 <pid>` on a running EEfinder now prints a full traceback of every
  thread (`faulthandler`), so a run that appears stuck can be diagnosed instead
  of guessed at.

### Notes
- The default `blastx` mode is the tested, reliable path. The DIAMOND modes can
  fail silently on newer `diamond` builds (subprocess stderr is discarded); verify
  the `diamond` build pinned in `env.yml` if a DIAMOND run produces no hits.

[2.0.0]: https://github.com/WallauBioinfo/EEfinder/compare/v1.1.2...v2.0.0
[1.1.2]: https://github.com/WallauBioinfo/EEfinder/releases/tag/v1.1.2
