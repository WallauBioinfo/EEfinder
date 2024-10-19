# EEfinder

EEfinder is a tool/python package that automatizes several tasks related to identification of Endogenous Elements present on Eukaryotic Genomes.

#### Install

EEfinder was developed and tested with BLAST 2.5.0 and DIAMOND 2.0.15, they are implemented on conda environments

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

#eefinder, version 1.1.1
```

For more information, check [EEfinder Wiki here](https://github.com/WallauBioinfo/EEfinder/wiki)

#### Cite us

If you use EEfinder in your research, please cite https://www.sciencedirect.com/science/article/pii/S2001037024003325:

```
Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). EEFinder: A general-purpose tool for identification of bacterial and viral endogenized elements in eukaryotic genomes. Computational and Structural Biotechnology Journal. https://doi.org/10.1016/j.csbj.2024.10.012
```