# Installation

EEfinder is a Python package that drives several external bioinformatics
binaries. Those binaries are **not** pip-installable, so the supported route is a
conda/mamba environment built from the bundled `env.yml`.

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

## Install with conda

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
# eefinder, version 1.1.1

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

## Next steps

- [Acquiring databases](databases.md) — the three reference inputs EEfinder
  needs.
- [Running EEfinder](running.md) — the example run against the bundled
  `test_files/`.
