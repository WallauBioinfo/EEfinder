# EEfinder

EEfinder is a tool/python package that automatizes several tasks related to identification of Endogenous Elements present on Eukaryotic Genomes.

## Endogenous Element Finder, Beta dev version

### Dependencies:

|Name|Version|
| ------ | ----- |
_libgcc_mutex|**0.1**
_openmp_mutex|**4.5**
bedtools|**2.27.1**
bio|**1.3.8**
biopython|**1.79**
biothings-client|**0.2.6**
blast|**2.5.0**
boost|**1.73.0**
bzip2|**1.0.8**
ca-certificates|**2022.3.29**
certifi|**2021.10.8**
charset-normalizer|**2.0.12**
icu|**58.2**
idna|**3.3**
ld_impl_linux-64|**2.35.1**
libboost|**1.73.0**
libffi|**3.3**
libgcc-ng|**9.3.0**
libgomp|**9.3.0**
libstdcxx-ng|**9.3.0**
lz4-c|**1.9.3**
mygene|**3.2.2**
ncurses|**6.3**
numpy|**1.22.3**
openssl|**1.1.1n**
pandas|**1.4.2**
pip|**21.2.4**
py-boost|**1.73.0**
python|**3.9.12**
python-dateutil|**2.8.2**
pytz|**2022.1**
readline**|8.1.2**
requests|**2.27.1**
setuptools|**61.2.0**
six|**1.16.0**
sqlite|**3.38.2**
tk|**8.6.11**
tqdm|**4.64.0**
tzdata|**2022a**
urllib3|**1.26.9**
wheel|**0.37.1**
xz|**5.2.5**
zlib|**1.2.11**
zstd|**1.4.9**

#### Install

EEfinder was developed and tested with BLAST 2.5.0 and DIAMOND 2.0.15, they are implemented on conda environments

```bash
git clone https://github.com/dezordi/EEfinder.git
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

### Running:

#### Test files:

There are some test files in the folder `test_data_dev`:

Argument|Test file
| --- | --- |
fasta_genome|Ae_aeg_Aag2_ctg_1913.fasta
virus_proteins_table|virus_subset2.csv
virus_db|virus_subset.fa
host_protein_db|filter_subset.fa

#### Test line:

```bash
eefinder -in test_files/Ae_aeg_Aag2_ctg_1913.fasta -mt test_files/virus_subset.csv -db test_files/virus_subset.fa -bt test_files/filter_subset.fa -od . -id
```
#### Default line:

```bash
eefinder -in <fasta_genome> -mt <virus_proteins_table> -db <virus_db> -bt <host_protein_db> -od . -id
```

#### Default line (DIAMOND):

You can choose which mode of DIAMOND you want to use between:
- fast
- mid-sensitive
- sensitive
- more-sensitive
- very-sensitive
- ultra-sensitive

Using 'fast' mode as example:

```bash
eefinder -in <fasta_genome> -mt <virus_proteins_table> -db <virus_db> -bt <host_protein_db> -od . -md fast -id
```

#### Keeping temporaries:

For not deleting temporaries after conclusion use the argument: `-rm`

#### Merge level:

For decide which philogetic level is going to be use, between family or genus, use the argument: `-ml --merge_level`

#### Modify min lenght of contigs:

In EEfinder you can choose the minimum lenght of contigs that BLAST/DIAMOND is going to use with the argument: `-ln --lenght`

#### Modify flank lenght:

You can choose the lenght of flanking regions going to be extracted with the argument: `-fl --flank`

#### Modify merge lenght:

There are the possibility to change the lenght for merging elements in EEfinder: `-lm --limit`

#### Name your prefix:

You can name the prefix that EEfinder is going to use to create your output files with: `-pr --prefix`
We suggest to use the name of the genome with the name of the assembly like: **Ae_aeg_Aag2 for Aedes _aegypti_/Aag2**

### Adquire datasets

#### Viral Protein database

For getting the **protein database** you can go on [ncbi virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and follow this steps:

Scroll down and select RefSeq
![Select RefSeq](images/protein_dataset/refseq.png)
Click on download button
![Click Download button](images/protein_dataset/download_button.png)
Select protein
![Select Protein](images/protein_dataset/select_protein.png)
Select Download all records
![Select Download all records](images/protein_dataset/download_all.png)
Select default
![Select default](images/protein_dataset/fasta_line.png)

#### Protein table

For getting the protein csv file continue on the RefSeq page and also click on the download button and follow this steps:
![Select csv](images/protein_table/csv.png)
Select csv
![Select Download all records](images/protein_table/download_all.png)
Select Download all records
![Follow the template](images/protein_table/template.png)

#### For bacterias

For bacterias you can use the addicional script **bac_retriever.py**, on the folder **accessory_scripts/**, that create the table with the proteins you choose.

```bash
python bac_retriever.py -in <fasta_with_bac_protein> -em <your_email_associated_to_ncbi_account> -key <your_ncbi_api_key>
```

This script will automatically recovery the accession ID present on fasta header and will create a tabular file with the same structure of the `virus_subset.tsv` test example, with the bacterial metadata information.

##### Disclaimer

Your bacterial protein fasta file **should** have sequences recovered from NCBI protein databases, and the sequences **should** have the same header pattern (starting with NCBI access protein code).

##### Email and APIkey

You have to register on NCBI and create an APIkey, for the APIkey go on **Account Settings** and find **API Key Management** to create your APIkey

![account_settings](images/apikey/account_settings.png)

## Filter Datasets

### Host proteins

The other dataset for cleaning is a set of host proteins, we suggest RefSeq proteins from NCBI.

## Outputs

Name|Meaning
| --- | --- |
prefix.EEs.fa|Fasta file with Endogenous Elements nucleotide sequences
prefix.EEs.tax.tsv|TSV file with Endogenous Elements taxonomy
prefix.EEs.flanks.fa|Fasta file with Endogenous Elements plus 10000nt in each flanking regions
prefix.EEs.cleaned.fa|Fasta file with Cleaned Endogenous Elements
prefix.EEs.cleaned.tax.tsv|TSV file with Cleaned Endogenous Elements

### Output diretory

```
outdir/
├── EEfinder.log.txt
├── prefix.EEs.cleaned.fa
├── prefix.EEs.cleaned.tax.tsv
├── prefix.EEs.fa
├── prefix.EEs.flanks.fa
├── prefix.EEs.tax.tsv
└── tmp_files
    ├── prefix.rn
    ├── prefix.rn.fmt
    ├── prefix.rn.fmt.blastx
    ├── prefix.rn.fmt.blastx.csv
    ├── prefix.rn.fmt.blastx.filtred
    ├── prefix.rn.fmt.blastx.filtred.bed
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.csv
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed
    ├── prefix.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank
    ├── prefix.rn.fmt.fai
    └── prefix.rn.fmt.rn.fmt.lenght
```
