# Overlap resolution (`--overlap`)

Two endogenous elements can occupy the same genomic region while being assigned
to **different families** — for example when a locus matches proteins from more
than one viral lineage. `TagElements` flags such elements in the
`tag` column (`overlaped`, versus `unique`), and `--overlap` (`-ov`) decides
what to do with them.

| Value | Behaviour |
|-------|-----------|
| `keep` (default) | Keep every element (no filtering). |
| `longest` | Among overlaping elements, keep the longest and drop the shorter ones. |
| `targets` | Resolve each overlap cluster by a family list — a keep-list (`--target_families`) or a drop-list (`--non_target_families`). |

## `targets`: keep-list vs drop-list

With `--overlap targets` you must provide **exactly one** of two repeatable
family lists:

- `--target_families` / `-tf` — **keep-list**: in a cluster that contains a
  target-family member, keep the target-family elements and drop the rest; a
  cluster with no target member is kept as-is.
- `--non_target_families` / `-ntf` — **drop-list**: drop the listed families from
  each cluster, unless *every* member of the cluster is a listed family (in which
  case the cluster is kept, never wiped).

```bash
# keep-list: keep only these families in an overlap
eefinder screening ... -ov targets -tf Flaviviridae -tf Caulimoviridae

# drop-list: remove these families from an overlap
eefinder screening ... -ov targets -ntf Retroviridae
```

Providing neither list (or both) with `--overlap targets` is an error.

## Where filtered elements go

Filtering is applied to **every** final result (`.EEs.fa`, `.EEs.tax.tsv`,
`.EEs.gff3`, `.EEs.flanks.fa`, and the `--clean_masked` variants). The
filtered-out elements are **not** discarded — they are written to a
`tmp_outputs/` directory (`PREFIX.EEs.removed.{fa,tax.tsv}`) for inspection, so
kept and removed elements together always partition the unfiltered run.

Overlap resolution runs **before** the GFF3 and flanking-region steps, so those
outputs only ever contain the kept elements.
