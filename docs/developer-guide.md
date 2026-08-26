# Developer guide

## Package layout

EEfinder uses a flat layout — the package lives in `eefinder/`, no `src/`. Each
processing step is a small **side-effect class** whose `__init__` runs the work
(files in, files out); there is no shared in-memory pipeline object. Steps
communicate through files on disk whose names accrete suffixes (`.rn.fmt`,
`.blastx`, `.filtred`, `.bed`, `.tax`, …).

| Module | Responsibility |
|--------|----------------|
| `scripts/main.py` | The `click` command group (`screening`, `get-databases`). |
| `prepare_data.py` | `PrepareGenome` — prefix FASTA headers **and** apply the length cutoff in one pass (what the pipeline runs); `InsertPrefix` standalone. |
| `clean_data.py` | `RemoveShortSequences`, `MaskClean`. |
| `make_database.py` | `MakeDB` — build BLAST/DIAMOND databases. |
| `similarity_analysis.py` | `SimilaritySearch` — the search (run twice). |
| `translation.py` | Protein prediction + coordinate traceback (`gv`/`rv`/`gv-rv`). |
| `filter_table.py` | `FilterTable` — redundant-hit filtering. |
| `bed.py` | `GetFasta`, `MergeBed`, `BedFlank`, … (bedtools wrappers). |
| `get_length.py` | `GetLength` — the `<id>\t<length>` index bedtools slop needs. |
| `compare_results.py` | `CompareResults` — host-bait filtering. |
| `get_taxonomy.py` | `GetTaxonomy` / `GetFinalTaxonomy` / `GetCleanedTaxonomy`. |
| `tag_elements.py` | `TagElements` — overlap flags + `Average_pident`. |
| `overlap.py` | `FilterOverlap` — overlap resolution strategies. |
| `gff.py` | `WriteGFF3`. |
| `get_databases.py` | `GetDatabases` — the `get-databases` implementation. |
| `taxon_exclusion.py` | Taxonomy expansion so `--exclude-taxon` can drop a branch before it is downloaded. |
| `progress.py` | Terminal progress reporting + download retry/stall detection. |
| `normalization.py` | `standardize_protein` — per-target protein-name cleaning. |
| `utils.py` | Path/timing helpers + the run-info dataclasses. |
| `versions.py` | Dependency-version detection + `env.yml` comparison. |
| `log.py` | The `eefinder` logger + `enable_debug()`. |
| `run_message.py` | The startup banner and citation notice. |

## Conventions & gotchas

- **Side-effect classes:** instantiating a step class runs it. Don't expect
  return values — check the output file.
- **Filename chaining:** downstream steps hard-code the accreted suffix of the
  upstream file. Changing an output name means updating every consumer in
  `main.py`.
- **The default `blastx` mode is the reliable path.** The DIAMOND modes can fail
  silently because the subprocess stderr is routed to `DEVNULL`; verify the
  `diamond` build (env pins `diamond=2.2.3`) if a DIAMOND run produces no hits.
- **Debug logging:** `--debug` (on both commands) lowers the `eefinder` logger to
  DEBUG via `log.enable_debug()`; the `logger.debug(...)` calls throughout are
  silent otherwise.

## Style

- Format with **black** (line length 88): `black eefinder tests`. The package and
  tests are black-clean; keep them that way.
- Target Python is pinned to **3.9** (`env.yml`). Avoid syntax that only parses
  on 3.10+ (e.g. parenthesised multi-item `with` statements); use
  `from __future__ import annotations` for newer typing syntax.
- Public functions/classes get NumPy-style docstrings.

## Build metadata

Build metadata lives entirely in `pyproject.toml`, built with the **hatchling**
backend. The `[project]` table declares the runtime dependencies, the `dev`
extra (`pytest` + `black`), and the `eefinder` console script; the same file
holds the black and pytest configuration. There is no `setup.py` or
`MANIFEST.in`. The version is set once in `[project].version` and read back at
runtime via stdlib `importlib.metadata` (no `pkg_resources` / `setuptools`
runtime dependency).

```bash
pip install .            # runtime install
pip install -e ".[dev]"  # editable install with pytest + black
python -m build          # build the wheel/sdist (needs the `build` package)
```

### Dependency bounds

Two upper bounds are deliberate and must not be relaxed without code changes:

| Pin | Reason |
|-----|--------|
| `biopython>=1.79,<1.86` | `Bio.Blast.Applications` (used by `make_database.py` and `similarity_analysis.py`) was **removed in Biopython 1.86**. Lifting the cap means replacing `NcbiblastxCommandline` / `NcbimakeblastdbCommandline` with plain `subprocess` calls. |
| `requires-python = ">=3.9,<3.12"` | Matches the interpreter range the pinned scientific stack supports. |

`pandas` and `numpy` are bounded only below the next major (`<3`): the unit
suite passes in full against pandas 2.3 / numpy 2.2 as well as the pandas 1.4 /
numpy 1.23 pinned in `env.yml`. `env.yml` pins the exact tested versions; the
`pyproject.toml` bounds describe the range the code is compatible with.

## Releasing to PyPI

Releases are published by `.github/workflows/publish.yml` when a GitHub release
is published. It uses **PyPI Trusted Publishing** (OIDC), so no API token is
stored in the repository; the workflow builds the sdist and wheel, runs
`twine check --strict`, installs the wheel into a clean virtualenv and runs
`eefinder --version` before uploading.

The one-time setup on PyPI (*Manage project → Publishing*, or *Add a pending
publisher* for the first ever release) registers the publisher as owner
`WallauBioinfo`, repository `EEfinder`, workflow `publish.yml`, environment
`pypi`.

To cut a release:

```bash
# 1. bump [project].version in pyproject.toml, commit
# 2. tag and publish a GitHub release -- the workflow does the rest
git tag v2.0.0 && git push origin v2.0.0
gh release create v2.0.0 --generate-notes
```

To check the artefacts by hand before trusting the workflow:

```bash
python -m build
twine check --strict dist/*
twine upload --repository testpypi dist/*     # TestPyPI first
```

```{warning}
A version published to PyPI is **immutable** — it cannot be overwritten or
re-uploaded, only yanked. Verify on TestPyPI before the real upload, and bump
the version rather than trying to replace a bad release.
```

## Continuous integration

`.github/workflows/tests.yml` runs the binary-free unit tests
(`pytest -m "not integration"`) on every pull request to `master`. It installs
only the pip runtime dependencies (kept in sync with `env.yml`) and pytest — no
BLAST/DIAMOND/bedtools — so the tool-gated tests skip cleanly. See
[Testing](testing.md) for the full suite.

## Building the docs locally

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
# open docs/_build/html/index.html
```

`docs/conf.py` reads the version straight out of `pyproject.toml`, so the docs
build does not need EEfinder — or BLAST, DIAMOND and bedtools — installed.

## Read the Docs

`.readthedocs.yaml` at the repository root drives the hosted build: Ubuntu
22.04, Python 3.11, Sphinx pointed at `docs/conf.py`, and the dependencies from
`docs/requirements.in`. Read the Docs builds one documentation version per git
branch or tag that it is configured to track, so the 1.x documentation stays
available from its own tag after later releases land on `master`.

`docs/requirements.txt` is a lock of the versions that resolve from
`requirements.in`; regenerate it with:

```bash
pip install -r docs/requirements.in && pip freeze > docs/requirements.txt
```
