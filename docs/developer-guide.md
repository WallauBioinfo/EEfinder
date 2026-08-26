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
| `utils.py` | `check_outdir`, `step_info`, `running_info` — path and run-log helpers. |
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
  [Acquiring databases](databases.md#the-metadata-csv-format). Adding a column to
  that file shifts every taxonomy assignment.
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

Build metadata lives in `setup.py` (setuptools). It declares the `eefinder`
console script (`eefinder.scripts.main:main`) and the `click` runtime pin; the
scientific stack (biopython, pandas, numpy) comes from `env.yml` rather than
from `install_requires`, because it must match the conda-provided binaries.

The version is set once in `setup.py` and read back at runtime by
`eefinder/__init__.py` through `pkg_resources`, which is why `eefinder --version`
reports the *installed* distribution version — reinstall after bumping it.

```bash
pip install .      # runtime install
pip install -e .   # editable install for development
```

## Building the docs locally

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
# open docs/_build/html/index.html
```

`docs/conf.py` reads the version straight out of `setup.py`, so the docs build
does not need EEfinder — or BLAST, DIAMOND and bedtools — installed.

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
