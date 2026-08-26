# EEfinder

EEfinder is a tool/python package that automatizes several tasks related to identification of Endogenous Elements present on Eukaryotic Genomes.

#### Install

EEfinder drives BLAST, DIAMOND and bedtools, which are not pip-installable — install them first, then the package.

##### From PyPI

```bash
conda create -n EEfinder -c conda-forge -c bioconda \
  "python>=3.9,<3.12" "blast>=2.5" "diamond>=2.0.15" "bedtools>=2.27"
conda activate EEfinder

pip install eefinder
```

##### From source

Cloning also gives you `env.yml`, which pins the exact versions EEfinder was developed and tested against (BLAST 2.5.0, DIAMOND 2.0.15, bedtools 2.27.1), plus the example data in `test_files/`.

```bash
git clone https://github.com/WallauBioinfo/EEfinder.git
cd EEfinder
conda env create -f env.yml
conda activate EEfinder
pip install .
```

#### Check tool

```bash
eefinder --version

#eefinder, version 1.1.2
```

#### Documentation

Full documentation is hosted on Read the Docs: **https://eefinder.readthedocs.io**

#### Cite us

If you use EEfinder in your research, please cite https://www.sciencedirect.com/science/article/pii/S2001037024003325:

```
Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). EEFinder: A general-purpose tool for identification of bacterial and viral endogenized elements in eukaryotic genomes. Computational and Structural Biotechnology Journal. https://doi.org/10.1016/j.csbj.2024.10.012
```