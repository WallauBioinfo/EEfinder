# Installation

EEfinder is a Python package that drives several external bioinformatics
binaries. Those binaries are **not** pip-installable, so a `pip install` alone
is never a complete installation — the conda/mamba environment built from the
bundled `env.yml` provides them.

## Requirements

| Dependency | Role | Provided by |
|------------|------|-------------|
| Python 3.9 | runtime | `env.yml` |
| BLAST 2.5.0 (`blastx`, `makeblastdb`) | similarity search + database build | `env.yml` (`blast`) |
| DIAMOND 2.0.15 (`diamond`) | fast alternative to BLAST | `env.yml` (`diamond`) |
| bedtools 2.27.1 | sequence extraction, merging, flank extraction | `env.yml` (`bedtools`) |
| biopython, pandas, numpy, click | Python runtime deps | `env.yml` (pip) |

EEfinder was developed and tested against the BLAST and DIAMOND versions pinned
above; `env.yml` installs exactly those.

## Install from PyPI

Install the external binaries first, then EEfinder itself:

```bash
# 1. the binaries (conda-forge / bioconda)
conda create -n EEfinder -c conda-forge -c bioconda \
  "python>=3.9,<3.12" "blast>=2.5" "diamond>=2.0.15" "bedtools>=2.27"
conda activate EEfinder

# 2. the package
pip install eefinder
```

## Install from source

Cloning gives you `env.yml` — which pins the exact versions EEfinder was tested
against — plus the example data in `test_files/`:

```bash
git clone https://github.com/WallauBioinfo/EEfinder.git
cd EEfinder

conda env create -f env.yml     # or: mamba env create -f env.yml
conda activate EEfinder

pip install .                   # or `pip install -e .` for development
```

## Verify the installation

```bash
eefinder --version
# eefinder, version 1.1.2

eefinder --help
# Usage: eefinder [OPTIONS]
#   This tool predict regions of Endogenous Elements in Eukaryote Genomes.
```

Then check that the external binaries are on `PATH` inside the activated
environment:

```bash
blastx -version
diamond --version
bedtools --version
```

```{note}
The Python dependencies are resolved by pip, but two of them are bounded on
purpose: **Biopython is capped below 1.86**, which removed the
`Bio.Blast.Applications` wrappers this version relies on, and Python is capped
below 3.12. Both bounds are declared in `pyproject.toml`, so pip will pick
compatible versions for you.
```

## Next steps

- [Acquiring databases](databases.md) — the three reference inputs EEfinder
  needs.
- [Running EEfinder](running.md) — the example run against the bundled
  `test_files/`.
