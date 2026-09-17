# Installation

EEfinder is a Python package that drives several external bioinformatics
binaries. Those binaries are **not** pip-installable, so a `pip install` alone
is never a complete installation — the binaries have to come from somewhere
else: the [Bioconda package](#install-from-bioconda), the
[container image](#run-with-a-container) (nothing to install at all), or a
conda/mamba environment (the bundled `env.yml`, or the one-liner under
[Install from PyPI](#install-from-pypi)).

## Requirements

| Dependency | Role | Provided by |
|------------|------|-------------|
| Python 3.9 | runtime | `env.yml`, Bioconda |
| BLAST 2.5.0 (`blastx`, `makeblastdb`) | similarity search + database build | `env.yml` (`blast`), Bioconda |
| DIAMOND 2.0.15 (`diamond`) | fast alternative to BLAST | `env.yml` (`diamond`), Bioconda |
| bedtools 2.27.1 | sequence extraction, merging, flank extraction | `env.yml` (`bedtools`), Bioconda |
| biopython, pandas, numpy, click | Python runtime deps | `env.yml` (pip), Bioconda, PyPI |

EEfinder was developed and tested against the BLAST and DIAMOND versions pinned
above; `env.yml` installs exactly those, and the Bioconda package requires at
least those versions. The container image carries all of them, including the
Python interpreter.

(install-from-bioconda)=
## Install from Bioconda

The [Bioconda package](https://anaconda.org/bioconda/eefinder) is the shortest
route: it declares BLAST, DIAMOND and bedtools as dependencies, so one command
installs EEfinder *and* the binaries it drives.

```bash
conda create -n EEfinder -c conda-forge -c bioconda eefinder
conda activate EEfinder
```

Or into an existing environment:

```bash
conda install -c conda-forge -c bioconda eefinder
```

```{tip}
`mamba` works as a drop-in replacement for `conda` in either command and
resolves considerably faster.
```

(run-with-a-container)=
## Run with a container

The Bioconda package is repackaged automatically as a
[BioContainers image](https://quay.io/repository/biocontainers/eefinder), which
bundles EEfinder with BLAST, DIAMOND and bedtools. Nothing is installed on the
host — useful on machines where you cannot create conda environments, and for
reproducible or HPC runs.

```{important}
BioContainers images are **always version-tagged and there is no `latest`
tag** — a tag has to be given explicitly. Check the
[tag list](https://quay.io/repository/biocontainers/eefinder?tab=tags) for the
current one; the examples below use `1.1.2--pyhdfd78af_0`, where the suffix is
the Bioconda build string.
```

### Docker

The container sees only what you mount into it, so the genome, the databases and
the output directory all have to live under a mounted path. Mounting the current
directory as the working directory is the simplest arrangement:

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$PWD":/data -w /data \
  quay.io/biocontainers/eefinder:1.1.2--pyhdfd78af_0 \
  eefinder --version
```

A full run from the repository root, over the bundled `test_files/` (the same
example as in [Running EEfinder](#example-run)):

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$PWD":/data -w /data \
  quay.io/biocontainers/eefinder:1.1.2--pyhdfd78af_0 \
  eefinder \
    -in test_files/Ae_aeg_Aag2_ctg_1913.fasta \
    -od results_test \
    -db test_files/virus_subset.fa \
    -mt test_files/virus_subset.csv \
    -bt test_files/filter_subset.fa \
    -ln 1000 -id -p 2 -lm 100
```

```{note}
`-u "$(id -u):$(id -g)"` makes the results belong to you instead of `root`. It
matters for more than tidiness: `-id` writes the BLAST/DIAMOND index files
*next to the database FASTA*, so that directory must be mounted writable by the
same user — index once outside the container, or keep the databases under the
mounted tree.
```

### Singularity / Apptainer

Pull the image once into a local `.sif` file:

```bash
singularity pull eefinder.sif \
  docker://quay.io/biocontainers/eefinder:1.1.2--pyhdfd78af_0
```

Or download the prebuilt Singularity image from the Galaxy depot, which skips
the Docker-to-SIF conversion:

```bash
curl -Lo eefinder.sif \
  https://depot.galaxyproject.org/singularity/eefinder:1.1.2--pyhdfd78af_0
```

Then run it. Singularity (and Apptainer, its current name — the commands are
interchangeable) bind-mounts the current directory and runs as you, so no `-v`
or `-u` equivalent is needed for data under `$PWD`:

```bash
singularity exec eefinder.sif eefinder --version

singularity exec eefinder.sif eefinder \
  -in test_files/Ae_aeg_Aag2_ctg_1913.fasta \
  -od results_test \
  -db test_files/virus_subset.fa \
  -mt test_files/virus_subset.csv \
  -bt test_files/filter_subset.fa \
  -ln 1000 -id -p 2 -lm 100
```

```{tip}
Paths outside the current directory need an explicit bind — databases on a
shared filesystem, for instance: `singularity exec -B /scratch/dbs:/dbs
eefinder.sif eefinder -db /dbs/viral.fa ...`.
```

(install-from-pypi)=
## Install from PyPI

The [PyPI package](https://pypi.org/project/eefinder/) ships the Python code
only. Install the external binaries first, then EEfinder itself:

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
# eefinder, version 1.1.3

eefinder --help
# Usage: eefinder [OPTIONS]
#   This tool predict regions of Endogenous Elements in Eukaryote Genomes.
```

With a container, prefix the same command — `docker run --rm
quay.io/biocontainers/eefinder:<tag> eefinder --version`, or `singularity exec
eefinder.sif eefinder --version` — and expect the version of the *image tag*,
which can lag the PyPI and Bioconda releases.

Then check that the external binaries are on `PATH` inside the activated
environment (inside a container they always are):

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
