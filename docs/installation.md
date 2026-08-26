# Installation

EEfinder is a Python package that drives several external bioinformatics
binaries. Those binaries are **not** pip-installable, so the supported route is a
conda/micromamba environment built from the bundled `env.yml`.

## Requirements

| Dependency | Role | Provided by |
|------------|------|-------------|
| Python 3.9 | runtime | `env.yml` |
| BLAST (`blastx`/`blastp`/`makeblastdb`) | similarity search + database build | `env.yml` (`blast`) |
| DIAMOND (`diamond`) | fast alternative to BLAST | `env.yml` (`diamond`) |
| bedtools | sequence extraction / merging | `env.yml` (`bedtools`) |
| NCBI datasets CLI (`datasets`) | database download (`get-databases`) | `env.yml` (`ncbi-datasets-cli`, pinned to 18.36.0) |
| cd-hit | dedup for `--translation_method gv-rv` **and** the `get-databases` `--cluster` step (on by default) | `env.yml` (`cd-hit`) |
| pyrodigal-gv / pyrodigal-rv | protein prediction (`gv`/`rv`/`gv-rv`) | `env.yml` (pip) |
| biopython (<1.86), pandas (<3), numpy (<3), click | Python runtime deps | `env.yml` (pip) |

## Install from PyPI

Install the external binaries first, then EEfinder itself:

```bash
# 1. the binaries (conda-forge / bioconda)
conda create -n EEfinder -c conda-forge -c bioconda \
  "python>=3.9,<3.12" "blast>=2.5" "diamond>=2.0.15" "bedtools>=2.27" \
  ncbi-datasets-cli cd-hit
conda activate EEfinder

# 2. the package
pip install eefinder
```

This is the quickest route, but it resolves the binaries to whatever versions
are current. To get the versions EEfinder is tested against, install from source
with `env.yml` as below.

## Install from source (recommended)

```bash
git clone https://github.com/WallauBioinfo/EEfinder.git
cd EEfinder

micromamba env create -f env.yml     # or: conda env create -f env.yml
micromamba activate EEfinder

pip install .                        # or `pip install -e .` for development
```

## Verify the installation

```bash
eefinder --version
# eefinder, version 2.0.0

eefinder --help
# Usage: eefinder [OPTIONS] COMMAND [ARGS]...
#   screening       Run the EEfinder screening pipeline on a genome.
#   get-databases   Download RefSeq protein databases (and metadata) ...
```

## Development install

To run the test suite and format the code, add the development dependencies:

```bash
pip install ".[dev]"                  # pytest + black
# or, equivalently:
pip install -r requirements-dev.txt
```

See the [Developer guide](developer-guide.md) and [Testing](testing.md) pages
for details.

```{tip}
When EEfinder is installed **outside** its source tree, the run log's
dependency-drift check needs to find the reference `env.yml`. Point it there
with `export EEFINDER_ENV_YML=/path/to/env.yml`.
```
