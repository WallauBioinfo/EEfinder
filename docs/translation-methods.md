# Translation methods (`--translation_method`)

By default the similarity search translates the genome in the six reading frames
with `blastx` / `diamond blastx`. `--translation_method` (`-tm`) swaps this for
**up-front protein prediction**, aligning the predicted proteins with `blastp` /
`diamond blastp`.

The chosen method applies to **both** searches in a run — the main EE search and
the host-bait search — because a single value is threaded into the
`SimilaritySearch` step. They therefore can never diverge.

## Methods

| Value | Behaviour |
|-------|-----------|
| `default` | Six-frame `blastx` / `diamond blastx` (current behaviour). |
| `gv` | Predict proteins with [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv), align with `blastp` / `diamond blastp`. |
| `rv` | Predict proteins with [pyrodigal-rv](https://github.com/LanderDC/pyrodigal-rv), align with `blastp` / `diamond blastp`. |
| `gv-rv` | Run **both** predictors, drop redundancy with `cd-hit` (100% identity / 100% coverage), then align. |

```bash
eefinder screening ... -tm gv
eefinder screening ... -tm gv-rv
```

## How the prediction path works

```{mermaid}
flowchart LR
    A[Nucleotide query] --> B[predict_proteins<br/>pyrodigal-gv / -rv]
    B --> C[Predicted-protein FASTA]
    B --> D[Coordinates TSV<br/>protein_id, contig, start, end, strand, tool]
    C --> E{gv-rv?}
    E -- yes --> F[cd-hit 100%/100%]
    E -- no --> G[blastp / diamond blastp]
    F --> G
    G --> H[traceback<br/>aa coords → nucleotide coords]
    D --> H
    H --> I["{query}.blastx<br/>same schema as blastx"]
```

Prediction (in `eefinder/translation.py`) writes, alongside the
predicted-protein FASTA, a **coordinates TSV**
(`protein_id, contig, start, end, strand, tool`). After the protein-vs-protein
search, a **coordinate traceback** maps each hit's amino-acid span back to
nucleotide coordinates on the source contig, following the `blastx` convention
(`qstart < qend` on the plus strand, `qstart > qend` on the minus strand).

Because the traceback emits the exact schema `blastx` would have produced,
`SimilaritySearch` always writes the same `{query}.blastx` table and **every
downstream step (filter / bed / taxonomy / GFF3) is unchanged** regardless of the
chosen method. The final Element-IDs are nucleotide coordinates in all modes.

```{note}
`gv-rv` uses `cd-hit` to drop proteins predicted identically by both tools.
Identical sequences have identical coordinates, so the coordinates TSV keeps
every `gv`+`rv` entry and the cluster representative id still resolves during
traceback.
```

## Dependencies

The prediction modes need extra dependencies, all pinned in `env.yml`:

- `pyrodigal-gv`, `pyrodigal-rv` (pip) — the predictors;
- `cd-hit` (conda) — only for `gv-rv`.

If they are absent, only `default` is available.

## Choosing a method

- **`default`** — the tested, reliable baseline; use it unless you have a reason
  not to.
- **`gv` / `rv`** — predict genes up front (viral-tuned Prodigal variants),
  which can behave better on gene-dense or well-assembled inputs. Note the gene
  boundaries pyrodigal calls will **not** line up 1:1 with six-frame `blastx`, so
  results are biologically different, not a drop-in reproduction.
- **`gv-rv`** — union of both predictors with redundancy removed; the broadest
  prediction-based recall.

```{warning}
The prediction modes change *which* proteins are searched, so their EE calls
differ from `default`. Validate the biological calls on a dataset you understand
before relying on them. The coordinate *mechanics* (traceback, both-search
consistency) are covered by the test suite and by a byte-for-byte golden
comparison for `gv-rv` (`test_files/expected_results/gv-rv/`).
```
