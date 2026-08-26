# EEfinder Documentation

**EEfinder** is a Python CLI and package for the identification of **Endogenous
Elements (EEs)** — sequences of non-host origin integrated into eukaryotic
genomes, most commonly derived from viruses or bacteria. A single command takes
a genome assembly and a reference protein database and returns the element
sequences, their taxonomic assignment, and their genomic flanking regions.

It was published in the *Computational and Structural Biotechnology Journal*
(Dias, Dezordi & Wallau, 2024,
[doi:10.1016/j.csbj.2024.10.012](https://doi.org/10.1016/j.csbj.2024.10.012)).

## Main features

Endogenous elements can be ancestral integrations, specially when we are looking
for Endogenous Virus Elements (EVEs). Once inserted, they are no
longer under selection to maintain viral or bacterial function, so they can
accumulate substitutions, indels and premature stop codons, and are frequently
truncated or split apart by subsequent genomic rearrangements. Detecting them is
involves a **remote-homology** problem, and each step of EEfinder addresses one
consequence of that.

- **Taxon-agnostic by design.** EEfinder encodes no donor taxon. It ships with
  no reference database, no taxon list and no clade-specific heuristics: *what*
  it searches for is set by the protein database (`-db`) and its metadata table
  (`-mt`), and what counts as host background is set by the bait set (`-bt`).
  The published applications are viral and bacterial endogenization, but the
  same pipeline screens for elements from any donor lineage whose proteins can
  be assembled into a FASTA — and, symmetrically, a search can be narrowed to a
  single family or genus simply by restricting the database. Every taxonomic
  name in the output is copied from the user's metadata table; the tool has no
  internal source of taxonomy. See
  [Acquiring databases](#databases-define-search).

- **Translated, protein-level similarity search.** Protein-based searches are generally
  more sensitive than nucleotide-based searches for detecting old viral or bacterial integrations,
  so EEfinder searches the genome translated in all six reading frames against a reference
  **protein** database. The search is parameterised for divergent sequences
  (BLOSUM45 substitution matrix, word size 3, *E*-value ≤ 1e-5), and alignments
  shorter than 33 aa are discarded. `blastx` is the default and the most
  sensitive option. `diamond blastx` is available as a faster alternative, but
  **recovers substantially fewer elements** — in the benchmark published with the tool,
  its most sensitive setting found 225 elements against BLAST's 481 on the same genome, with the
  loss concentrated on the most divergent sequences. See
  [The sensitivity trade-off](#diamond-sensitivity).

- **Host-gene bait filtering.** Reference viral and bacterial proteins share
  conserved domains with host proteins — reverse transcriptases, helicases,
  polymerases — and many host genes are themselves of retroelement origin. A
  locus can therefore hit a viral protein without being virus-derived. Every
  candidate is re-searched against a user-supplied set of host proteins (`-bt`)
  and retained only if its hit to the viral/bacterial database outscores its
  best host hit.

- **Redundant-hit resolution.** A single locus typically yields many overlapping
  alignments against many homologous references. Hits are collapsed per query,
  coordinate range and strand to the single highest-scoring alignment
  (`--range_junction`), so element counts reflect genomic loci rather than
  database redundancy.

- **Reconstruction of fragmented elements.** Deletions and rearrangements mean
  one ancestral integration often survives as several fragments, each with a
  slightly different — or truncated — taxonomic assignment. Neighbouring
  fragments on the same strand that share a taxon at the chosen rank are merged
  into one element (`--limit`, `--merge_level` at `genus` or `family`), and the
  contributing accessions, products and taxa are carried into the output.

- **Repeat-aware filtering.** EEs are enriched in transposable-element-rich,
  repetitive genomic compartments, which assemblies commonly deliver
  soft-masked. `--clean_masked` writes an additional set of outputs restricted
  to elements whose soft-masked and ambiguous base content stays below
  `--mask_per`, alongside the unfiltered set.

- **Explicit treatment of ambiguous assignments.** Because conserved domains are
  shared across viral families, one locus can be supported by references from
  more than one family. EEfinder does not silently pick a winner: overlapping
  elements assigned to different families are cross-referenced and flagged
  (`overlaped` / `unique`), and every element carries its supporting accessions
  with their individual percent identities, so
  the evidence behind an assignment can be inspected.

- **Insertion context.** `--flank` nt of genomic sequence are extracted on each
  side of every element, for characterising the integration site — neighbouring
  transposable elements, host genes and the integration junction itself.

- **Auditable runs.** A JSON `eefinder.log` records the resolved arguments, the
  parameters used at each step and per-step timing; by default every
  intermediate file is preserved under `tmp_files/`, so any step of the analysis
  can be re-examined.

## Documentation contents

```{toctree}
:maxdepth: 2

installation
databases
running
custom-arguments
output
accessory-scripts
developer-guide
```

## Quick links

- [GitHub repository](https://github.com/WallauBioinfo/EEfinder)
- [Issues / bugs](https://github.com/WallauBioinfo/EEfinder/issues)

## Cite us

If you use EEfinder in your research, please cite:

> Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). EEfinder: A
> general-purpose tool for identification of bacterial and viral endogenized
> elements in eukaryotic genomes. *Computational and Structural Biotechnology
> Journal*. https://doi.org/10.1016/j.csbj.2024.10.012
