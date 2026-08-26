# EEfinder bundled data

## `viral_proteins.tsv` — viral protein-name standardization map

A starting reference to normalise the free-text `Protein` names found in the
viral database metadata (built by `eefinder get-databases virus`) into a small
set of canonical names.

This is intentionally **not** exhaustive — it covers the most frequent /
biologically relevant proteins and is meant to be extended.

The map is loaded by `eefinder/normalization.py` (`standardize_protein`, target
`virus`).

### Columns

```
suggested_name <TAB> current_name <TAB> match_type <TAB> molecule_type_scope <TAB> notes
```

- **suggested_name** — the canonical name emitted when the rule matches.
- **current_name** — the already-normalised key to match against (lower-case,
  separators collapsed; see below).
- **match_type** — `exact` or `contains` (see step 2).
- **molecule_type_scope** — restricts the rule to a genome group (see shorthand).
- **notes** — free-text documentation for the rule (optional).

### How the table is applied

1. **Normalise** the raw `Protein` string before matching:
   - remove leaked NCBI `[key=value]` tags (e.g. `[organism=...]`);
   - strip a leading `CDS:` / `ORF:` naming directive;
   - normalise molecular-weight tokens (`33 kDa`, `33-kDa`, `33K-like protein`
     → `33 kDa protein`);
   - fix common misspellings (e.g. `membran` → `membrane`);
   - lowercase;
   - replace hyphens/underscores AND compound separators `/\();,` with a single
     space (so `CP/RdRp fusion`, `...; RdRp` become matchable tokens);
   - collapse runs of whitespace to one space; strip ends;
   - strip leading qualifiers: `putative`, `predicted`, `probable`, `possible`,
     `presumed`, `presumptive` (also from the emitted name);
   - strip trailing `, partial` / ` precursor`.

   e.g. `Putative RNA-dependent RNA Polymerase` → `rna dependent rna polymerase`.

2. **Match** against `current_name` (already normalised in this file):
   - `match_type=exact` — the normalised string equals `current_name`;
   - `match_type=contains` — `current_name` occurs as a whole-word substring.

   Prefer `exact` for short/ambiguous names (e.g. `l`, `n`, `cp`) so that, for
   example, `rna polymerase sigma factor` (a phage enzyme) does **not** become
   RdRp. On the **output** name, quotes and the characters `:,/\?!` are removed,
   the first letter is capitalised, and a name that is only a directive becomes
   `Unknown`.

3. **Scope** — only apply the row when the record's `Molecule_type` is within
   `molecule_type_scope`; rows scoped to a group must not rewrite proteins from
   other groups.

### `molecule_type_scope` shorthand

| Token    | Meaning                                                        |
|----------|----------------------------------------------------------------|
| `RNA`    | any RNA genome, RT excluded [ssRNA(+), ssRNA(-), ssRNA(+/-), dsRNA] |
| `+ssRNA` | ssRNA(+)                                                       |
| `-ssRNA` | ssRNA(-)                                                       |
| `dsRNA`  | dsRNA                                                          |
| `RT`     | reverse-transcribing [ssRNA-RT, dsDNA-RT]                     |
| `dsDNA`  | dsDNA                                                          |
| `ssDNA`  | ssDNA                                                          |
| `any`    | no restriction                                                |

### Row groupings (order in the file)

- **RNA viruses** — RdRp (the worked example: the viral replicase / L protein;
  `contains` catches compound names like `P2-RdRp`, `RdRp protein`, `CP/RdRp
  fusion`, `...; RdRp`, `RdRp-like`), Nucleocapsid Protein, Phosphoprotein (P),
  Matrix protein (M), Glycoprotein (envelope glycoprotein; covers RNA **and**
  DNA viruses), Fusion (F) — kept **separate** from the generic glycoprotein,
  Capsid / coat protein (CP), and other (+)RNA proteins.
- **RT viruses** — Reverse Transcriptase, Gag, Env, integrase.
- **dsDNA phages / large DNA viruses** — dominant in a whole-virome download;
  structural naming is a large separate domain, so this is a representative
  starter set to be expanded. Note: in dsDNA phages the major capsid protein
  (MCP) is a different context from the (+)RNA/dsRNA capsid protein.
- **Unresolved** — leading `CDS:`/`ORF:` naming directives are stripped in code
  before matching; a name that is only such a directive (e.g. `CDS:`, `ORF`)
  becomes `Unknown`. Names beginning with `hypothetical` (any spelling) are also
  flagged `Unknown` and dropped from the FASTA and CSV.

## `ictv_genome_composition.tsv`

ICTV genome-composition table (family → molecule type), sourced from
<https://ictv.global/virus-properties>. Used to fill the `Molecule_type` column
of the metadata CSV, which the NCBI datasets report does not provide.
