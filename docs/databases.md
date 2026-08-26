# Acquiring databases

A run needs three reference inputs:

| Input | CLI flag | Where it comes from |
|-------|----------|---------------------|
| Protein database FASTA | `-db` | NCBI Virus (viruses) or your own selection (bacteria) |
| Protein metadata CSV | `-mt` | NCBI Virus download, or `bac_retriever.py` |
| Host-gene bait FASTA | `-bt` | RefSeq proteins of the host (or a related) species |

(databases-define-search)=
## The databases define the search space

These three files are not just inputs — they are what makes the run a search for
one thing rather than another. EEfinder has no built-in reference data and no
notion of which taxa are "interesting": it reports similarity to whatever proteins
you give it, labelled with whatever taxonomy your metadata table carries.

Three consequences are worth stating explicitly:

- **The donor lineage is your choice.** A RefSeq viral protein set screens for
  endogenous viral elements; a bacterial protein set screens for endogenous
  bacterial elements. Nothing in the pipeline is specific to either, so any
  donor lineage whose proteins you can assemble into a FASTA can be screened
  for on the same code path.
- **Scope follows the database.** Restricting `-db` to a single family — or a
  single genus — restricts the whole analysis to it, at a fraction of the
  runtime. Conversely, a hit can only ever be reported for a lineage that is
  represented in the database: absence of an element is evidence of absence
  *from your reference set*, not from the genome.
- **The taxonomy is yours too.** `Family`, `Genus`, `Species`, `Molecule_type`
  and `Host` in the output tables are copied verbatim from `-mt`. EEfinder never
  consults an external taxonomy, so the rank names, spellings and any
  `Unclassified` entries in your results are the ones in your table.

The `-bt` bait set defines the opposite side of the comparison: it is the
background that candidate elements must *fail* to match better than they match
`-db`. It should therefore represent the host's own proteome (see
[Host-gene baits](#host-gene-baits)).

```{note}
The one taxonomic assumption the pipeline does make is about **rank**, not
lineage: `--merge_level` merges fragmented elements at `genus` or `family`, so
those two columns must be populated in `-mt` for merging to work as intended.
```

## The metadata CSV format

`-mt` is the table that turns a protein accession into a taxonomic assignment.
EEfinder reads it **by column position**, so it must have exactly these seven
columns, in this order, with a header row:

```text
Accession,Species,Genus,Family,Molecule_type,Protein,Host
YP_009664712.1,Bas-Congo tibrovirus,Tibrovirus,Rhabdoviridae,ssRNA(-),N protein,Homo sapiens
YP_009665181.1,Chick syncytial virus,Gammaretrovirus,Retroviridae,ssRNA-RT,polymerase,
```

The `Accession` values must match the first token of the corresponding FASTA
headers in `-db` (e.g. `>YP_009664712.1 N protein [Bas-Congo tibrovirus]`).
`test_files/virus_subset.csv` is a working miniature example.

```{warning}
Extra, missing or reordered columns will not raise an error — they will produce
a silently wrong taxonomy table, because the columns are addressed by index.
Check the header before a full run.
```

## Viral database and metadata

Both files come from the same [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/)
page, downloaded twice: once as FASTA (the database) and once as CSV (the
metadata).

### Protein database (`-db`)

Scroll down and select **RefSeq protein**:

![NCBI Virus RefSeq protein](_static/images/db_refseq.png)

Click the **Download** button:

![Download button](_static/images/db_download_button.png)

Select **Protein** as the sequence type:

![Select protein](_static/images/db_select_protein.png)

Select **Download all records**:

![Download all records](_static/images/db_download_all.png)

Keep the **default** FASTA definition line:

![Default FASTA definition line](_static/images/db_fasta_line.png)

### Protein table (`-mt`)

Staying on the same RefSeq page, click **Download** again and this time choose
**CSV format**:

![Select CSV](_static/images/table_csv.png)

Select **Download all records**:

![Download all records](_static/images/table_download_all.png)

Then select exactly the columns EEfinder expects — **Accession** (*with
version*), **Species**, **Genus**, **Family**, **Molecule type**, **Protein**
and **Host** — and nothing else:

![Column selection](_static/images/table_template.png)

## Bacterial database and metadata

Assemble the protein FASTA yourself from the NCBI protein database, then build
the matching metadata table with the accessory script
[`accessory_scripts/bac_retriever.py`](https://github.com/WallauBioinfo/EEfinder/blob/master/accessory_scripts/bac_retriever.py):

```bash
python accessory_scripts/bac_retriever.py \
  -in <fasta_with_bacterial_proteins> \
  -em <your_ncbi_email> \
  -key <your_ncbi_api_key>
```

The script reads the accession IDs from the FASTA headers and queries Entrez for
the metadata, writing a table with the same seven-column structure as the viral
example.

```{important}
The bacterial protein FASTA **must** come from an NCBI protein database, and its
headers **must** start with the NCBI protein accession — that is how the script
finds the records.
```

### Email and API key

Entrez requires an email address, and an API key raises the request rate limit.
Register at NCBI, then open **Account Settings** and find **API Key Management**:

![NCBI API key management](_static/images/ncbi_apikey.png)

(host-gene-baits)=
## Host-gene baits (`-bt`)

The bait set is a collection of **host proteins** used to discard false
positives: a putative EE that matches a host gene with a higher bitscore than it
matches the viral/bacterial reference is removed (see
[Running EEfinder](running.md#pipeline-steps)).

We suggest the RefSeq proteins of the host species, or of the closest available
relative. No metadata table is needed for this input — only the FASTA.
`test_files/filter_subset.fa` shows the expected shape (RefSeq *Aedes
albopictus* proteins used against an *Aedes aegypti* genome).

## Indexing

Whichever route you took, the databases must be indexed once for the search
engine you intend to use. Pass `-id`/`--index_databases` on the first run against
a given pair of databases and EEfinder will build the BLAST (`makeblastdb`) or
DIAMOND (`diamond makedb`) index for both `-db` and `-bt`. Later runs against the
same files can omit the flag.
