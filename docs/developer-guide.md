# Developer guide

## Package layout

EEfinder uses a flat layout — the package lives in `eefinder/`, no `src/`. Each
processing step is a small **side-effect class** whose `__init__` runs the work
(files in, files out); there is no shared in-memory pipeline object. Steps
communicate through files on disk whose names accrete suffixes (`.rn`, `.fmt`,
`.blastx`, `.filtred`, `.bed`, `.tax`, …).

| Module | Responsibility |
|--------|----------------|
| `scripts/main.py` | The `click` command, its options and the step orchestration. |
| `prepare_data.py` | `InsertPrefix` — prefix the FASTA headers with `PREFIX/`. |
| `clean_data.py` | `RemoveShortSequences` (the `--length` cutoff), `MaskClean` (the `--mask_per` filter). |
| `make_database.py` | `MakeDB` — build BLAST or DIAMOND databases. |
| `similarity_analysis.py` | `SimilaritySearch` — the search (run twice: database, then host baits). |
| `filter_table.py` | `FilterTable` — redundant-hit filtering and sense assignment. |
| `bed.py` | `GetFasta`, `GetAnnotBed`, `MergeBed`, `RemoveAnnotation`, `BedFlank`, `GetBed` (bedtools wrappers). |
| `get_length.py` | `GetLength` — the `<id>\t<length>` index `bedtools slop` needs. |
| `compare_results.py` | `CompareResults` — host-bait filtering. |
| `get_taxonomy.py` | `GetTaxonomy` / `GetFinalTaxonomy` / `GetCleanedTaxonomy`. |
| `tag_elements.py` | `TagElements` — overlap flags + `Average_pident`. |
| `utils.py` | `check_metadata_columns`, `check_metadata_file`, `check_outdir`, `step_info`, `running_info` — input-validation, path and run-log helpers. |
| `log.py` | The `eefinder` logger. |
| `run_message.py` | The startup banner and the citation notice. |

## Conventions & gotchas

- **Side-effect classes:** instantiating a step class runs it. Don't expect
  return values — check the output file.
- **Filename chaining:** downstream steps hard-code the accreted suffix of the
  upstream file. Changing an output name means updating every consumer in
  `main.py`.
- **The metadata CSV is read by column position.** `GetFinalTaxonomy` indexes
  into the joined table at fixed offsets (15–20), which assumes the seven-column
  `-mt` layout documented in
  [Acquiring databases](databases.md#the-metadata-csv-format). That layout is
  enforced rather than assumed: `check_metadata_columns` in `utils.py` holds the
  canonical column list, errors on a missing column, and warns and reorders (or
  drops extras) so `GetTaxonomy` always joins the columns in the order those
  offsets expect. `check_metadata_file` runs the same check on the header alone,
  from `main.py`, so a bad `-mt` fails before the similarity search. Changing the
  expected columns means changing that list *and* the offsets together.
- **The default `blastx` mode is the reliable path.** The DIAMOND modes can fail
  silently because the subprocess stderr is routed to `DEVNULL`; verify the
  `diamond` build (env pins `diamond=2.0.15`) if a DIAMOND run produces no hits.
- **`FilterTable` writes chunks to `{outdir}/tmp/`** while it works and removes
  that directory afterwards, so two runs must not share an output directory
  concurrently.

## Style

- Format with **black** (default line length): `black eefinder`.
- Target Python is pinned to **3.9** (`env.yml`). Avoid syntax that only parses
  on 3.10+ (e.g. parenthesised multi-item `with` statements).
- Public functions/classes carry a short docstring listing their keyword
  arguments.

## Build metadata

Build metadata lives entirely in `pyproject.toml`, built with the **hatchling**
backend. The `[project]` table declares the runtime dependencies, the `eefinder`
console script and the project URLs. There is no `setup.py` and no
`MANIFEST.in`. The version is set once in `[project].version` and read back at
runtime through stdlib `importlib.metadata`, so `eefinder --version` reports the
*installed* distribution version — reinstall after bumping it.

```bash
pip install .            # runtime install
pip install -e .         # editable install for development
python -m build          # build the wheel/sdist (needs the `build` package)
```

### Dependency bounds

Two upper bounds are deliberate and must not be relaxed without code changes:

| Pin | Reason |
|-----|--------|
| `biopython>=1.79,<1.86` | `Bio.Blast.Applications` (used by `make_database.py` and `similarity_analysis.py`) was **removed in Biopython 1.86**. Lifting the cap means replacing `NcbiblastxCommandline` / `NcbimakeblastdbCommandline` with plain `subprocess` calls. |
| `requires-python = ">=3.9,<3.12"` | Matches the interpreter range the pinned scientific stack supports. |

`pandas` and `numpy` are bounded only below the next major (`<3`); the
data-processing modules were verified against pandas 2.3 and numpy 2.4.

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
git tag v1.1.3 && git push origin v1.1.3
gh release create v1.1.3 --generate-notes
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

## The Bioconda package

EEfinder is also distributed on Bioconda
([`bioconda::eefinder`](https://anaconda.org/bioconda/eefinder)), whose recipe
lives in
[`bioconda-recipes/recipes/eefinder`](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/eefinder)
— not in this repository. It is the install route the documentation
recommends, because on top of the Python dependencies its `run:` section
requires `blast`, `diamond` and `bedtools`, so the binaries come with the
package.

Once a recipe build lands, Bioconda also publishes the matching
[BioContainers image](https://quay.io/repository/biocontainers/eefinder)
automatically — one tag per version *and build number*
(`<version>--<build string>`), never a `latest` — so the container route needs
no work in this repository either.

A new PyPI release is normally picked up by Bioconda's autobump bot, which opens
a PR bumping the version and the sdist hash. Two things still need attention
there:

- **Dependency changes do not propagate.** If a release changes the bounds in
  `[project].dependencies` or `requires-python`, the recipe's `run:` section has
  to be edited to match, in a PR to `bioconda-recipes`.
- **The version there lags.** Until that PR is merged, `conda install` still
  serves the previous release while PyPI already has the new one.

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
