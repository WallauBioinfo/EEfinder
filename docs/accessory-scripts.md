# Accessory scripts

[`accessory_scripts/`](https://github.com/WallauBioinfo/EEfinder/tree/master/accessory_scripts)
holds standalone helpers that surround the pipeline: preparing its inputs, and
summarising or plotting its outputs. They are **not** installed by
`pip install .` — run them from the cloned repository with the `EEfinder`
environment active.

```{note}
These scripts are companions rather than part of the supported pipeline. Several
of them shell out to tools that are not in `env.yml` (`cd-hit-est`, `mafft`,
`getorf` from EMBOSS), and the older ones still refer to EEfinder's former name,
*PEVEI*, in their help text. Install the extra tools yourself if you need them.
```

## Preparing inputs

### `bac_retriever.py`

Builds the `-mt` metadata table for a bacterial protein FASTA by querying Entrez
for each accession found in the headers. See
[Acquiring databases](databases.md#bacterial-database-and-metadata).

```bash
python accessory_scripts/bac_retriever.py -in <proteins.fa> -em <email> -key <api_key>
```

### `ictv_ncbi.py`

Takes a three-column CSV and returns protein information grouped by viral
family, for assembling family-scoped protein sets.

## Post-processing elements

### `get_copies.py`

Clusters the elements in a FASTA with `cd-hit-est` and reformats the cluster
output, so repeated copies of the same integration can be counted.

```bash
python accessory_scripts/get_copies.py -in <PREFIX.EEs.fa> -p 4 -m 4000
```

### `get_nr_orfs.py`

Runs EMBOSS `getorf` over a FASTA and removes redundant ORFs by clustering the
result — useful for inspecting the coding potential retained by an element.

```bash
python accessory_scripts/get_nr_orfs.py -in <PREFIX.EEs.fa> -ms 300 -p 4 -m 4000
```

### `CIAling_auto.py`

Clusters a FASTA with `cd-hit-est` and formats the output for downstream
alignment (`mafft`).

### `all_tax_to_tsv.py`

Takes a list file naming several `*.EEs.tax.tsv` files (one per assembly) and
concatenates them into a single TSV organised by assembly — the input the
plotting scripts below expect.

```bash
python accessory_scripts/all_tax_to_tsv.py -in tax.lst
```

## Plotting

| Script | Produces |
|--------|----------|
| `tax_to_barplot.py` | Stacked barplots of taxon composition, one per assembly, from a concatenated taxonomy file. |
| `tax_to_beeswarm.py` | Beeswarm plots of the elements in a `.tax` file, by family/order and protein category. |
| `flank_TEs.py` | Stacked barplots relating elements to the transposable elements found in their flanking regions. |
| `eefinder_plots/taxonomy_plots.R` | R/ggplot2 summaries (molecule type, family and genus counts) of a single `PREFIX.EEs.tax.tsv`. |

The Python plotting scripts share a common shape: `-in` for the input table,
`-md`/`--mode` to choose the taxonomic level (`Family` by default, or `Order`),
`-st`/`--specifictax` to restrict the plot to named taxa, and `-gt`/`--groupbytax`
to group rare taxa into an *others* category (default: fewer than 5 elements).

```bash
python accessory_scripts/tax_to_barplot.py -in all_taxonomy.tsv -md Family -gt 5
python accessory_scripts/tax_to_beeswarm.py -in PREFIX.EEs.tax.tsv -st Rhabdoviridae Flaviviridae
```

`taxonomy_plots.R` reads its input filename from the top of the script — edit
that line to point at your own taxonomy table before running it.
