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

#eefinder, version 1.0.0
```

For more information, check [EEfinder Wiki here](https://github.com/WallauBioinfo/EEfinder/wiki)
