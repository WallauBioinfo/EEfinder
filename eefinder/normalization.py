"""Protein-name standardisation used by ``get-databases``.

Standardisation collapses the free-text protein products in a RefSeq download to
a small set of canonical names.  Every target type shares a generic cleaning
pipeline (directive stripping, molecular-weight and misspelling normalisation,
special-character removal, capitalisation, and the bare-``CDS``/``ORF`` ->
``"Unknown"`` collapse); a target-specific *mapper* is layered on top:

* ``virus`` — matches the bundled ``data/viral_proteins.tsv`` map, respecting the
  molecule-type scope (e.g. every RdRp spelling/synonym -> ``RdRp``).
* ``bacteria`` — generic cleaning only for now (extension point for a future
  bacterial protein map).
* ``host`` — generic cleaning only (host baits are gene/protein names kept
  as-is aside from cleaning).

:func:`standardize_protein` is the public entry point; it dispatches on the
``target`` argument.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, Optional

# ---------------------------------------------------------------------------
# Shared cleaning primitives (target-agnostic)
# ---------------------------------------------------------------------------

#: Regex to capture molecular weight-like names such as "100 kDa", "33-kDa",
#: "33K-like protein", "33L protein" and standardize them to "X kDa protein".
#: The unit may be separated from the number by whitespace and/or a hyphen so
#: that "33 kDa" and "33-kDa" collapse to the same "33 kDa protein".
_MW_RE = re.compile(
    r"\b(\d+)[\s-]*(?:kda|kd|k|l)(?:[- ]*(?:putative\s+)?(?:nonstructural\s+)?(?:protein|like\s+protein))?\b",
    re.IGNORECASE,
)


def normalize_molecular_weight(name: str) -> str:
    """Standardise molecular-weight protein names to ``"X kDa protein"``.

    Parameters
    ----------
    name : str
        Raw protein product possibly carrying a molecular-weight token.

    Returns
    -------
    str
        ``name`` with any molecular-weight token rewritten to ``"X kDa
        protein"``; unchanged when no such token is present.

    Examples
    --------
    >>> normalize_molecular_weight("100 kDa")
    '100 kDa protein'
    >>> normalize_molecular_weight("33-kDa")
    '33 kDa protein'
    >>> normalize_molecular_weight("33K-like protein")
    '33 kDa protein'
    >>> normalize_molecular_weight("33KD putative nonstructural protein")
    '33 kDa protein'
    """
    return _MW_RE.sub(r"\1 kDa protein", name)


#: NCBI CDS FASTA "[key=value]" metadata tags (e.g. "[organism=...]",
#: "[protein=...]", "[gbkey=CDS]") that can leak into a raw product name. The
#: value may itself contain one level of nested "[...]" (a strain/isolate tag,
#: e.g. "[organism=Maize streak virus - A[South Africa]]").
_BRACKET_TAG_RE = re.compile(r"\[[A-Za-z_]+=(?:[^\[\]]|\[[^\[\]]*\])*\]")


def strip_bracket_tags(name: str) -> str:
    """Remove NCBI ``[key=value]`` metadata tags from a protein name.

    Parameters
    ----------
    name : str
        A protein product possibly carrying leaked ``[key=value]`` tags such as
        ``[organism=...]``, ``[protein=...]`` or ``[gbkey=CDS]``.

    Returns
    -------
    str
        ``name`` with any such tag removed (whitespace is not collapsed here).

    Examples
    --------
    >>> strip_bracket_tags("nucleoprotein [organism=Rabies lyssavirus]").strip()
    'nucleoprotein'
    """
    return _BRACKET_TAG_RE.sub(" ", name)


#: Special characters (incl. quotes) stripped from a protein name.
_SPECIAL_CHARS = ":,/\\?!\"'"

_QUALIFIER_RE = re.compile(r"^(putative|predicted|probable|hypothetical)\s+")
_TRAILING_RE = re.compile(r"(,\s*partial|\s+precursor)$")

#: Leading hedging qualifiers removed from the *emitted* name (case-insensitive,
#: repeated). "hypothetical" is deliberately excluded so "hypothetical protein"
#: does not collapse to a bare "protein".
_OUTPUT_QUALIFIER_RE = re.compile(
    r"^\s*(?:putative|putatively|predicted|probable|possible|presumed|presumptive)\s+",
    re.IGNORECASE,
)


def _strip_qualifiers(text: str) -> str:
    """Drop leading hedging qualifiers (e.g. "putative ") from a protein name."""
    prev = None
    while prev != text:
        prev = text
        text = _OUTPUT_QUALIFIER_RE.sub("", text)
    return text


#: Non-structural-protein designation (``NS5``, ``NS5A``, ``NS4B``, ...) followed
#: by a redundant "protein" / "-like protein" / "peptide" suffix.  The bare
#: designation is kept (e.g. ``NS5 protein`` / ``NS5-like protein`` -> ``NS5``).
_NS_DESIGNATION_RE = re.compile(
    r"^(NS\d+[A-Za-z]?)(?:[\s-]+(?:like[\s-]+)?protein|[\s-]+peptide)$",
    re.IGNORECASE,
)


def _strip_designation_suffix(text: str) -> str:
    """Drop a redundant "protein"/"peptide" suffix from an ``NSxx`` designation."""
    match = _NS_DESIGNATION_RE.match(text.strip())
    return match.group(1).upper() if match else text


#: Leading "CDS:" / "ORF:" naming directives to strip from a protein name.
_LEADING_DIRECTIVE_RE = re.compile(r"^\s*(cds|orf)\s*:\s*", re.IGNORECASE)

#: Names that carry no protein information once directives/punctuation are gone.
_UNKNOWN_TOKENS = {"", "cds", "orf"}

#: Common misspellings/truncations in NCBI protein names, observed in a full
#: RefSeq viral download.  Applied by :func:`_apply_typos` both to the match key
#: (so the corrected form matches a protein map) and to the emitted name (so typo
#: variants of an *unmapped* protein still converge).  Keys are matched as whole
#: words, case-insensitively (see :data:`_TYPO_RE`), so a typo that is a prefix
#: of the correct spelling (e.g. "membran" -> "membrane") does not corrupt the
#: already-correct word.
_TYPO_CORRECTIONS: dict[str, str] = {
    # Truncations/mangled "polymerase".  "polymeras" only matches when it is not
    # already the correct "polymerase" (word boundary requires no trailing "e").
    "polymeras": "polymerase",
    "polymrease": "polymerase",
    "polymarase": "polymerase",
    "polymerse": "polymerase",
    "polymease": "polymerase",
    "polyermase": "polymerase",
    # Nucleocapsid / capsid.
    "nucleocapside": "nucleocapsid",
    "nucleopasid": "nucleocapsid",
    "capside": "capsid",
    "caspsid": "capsid",
    # Polyprotein.
    "polyprotien": "polyprotein",
    "polyportein": "polyprotein",
    "plyprotein": "polyprotein",
    # Membrane / phospho.
    "membran": "membrane",
    "membraine": "membrane",
    "membrain": "membrane",
    "phoshoprotein": "phosphoprotein",
    # Hypothetical (so the "hypothetical*" drop below catches misspellings).
    "hypotheticla": "hypothetical",
    "hypotheticl": "hypothetical",
    "hyphothetical": "hypothetical",
    "hypotetical": "hypothetical",
    "hypothecial": "hypothetical",
    "hyppothetical": "hypothetical",
    # Truncated "glycoprotein" (e.g. "Glycop C", "Glycoprot").  Word boundaries
    # keep the already-correct "glycoprotein"/"glycoproteins" untouched.
    "glycoprot": "glycoprotein",
    "glycop": "glycoprotein",
    # Compound words where "glycoprotein" is fused with a prefix — split so the
    # ``\bglycoprotein\b`` contains-match can find them.
    "phosphoglycoprotein": "phospho glycoprotein",
    "proteinglycoprotein": "protein glycoprotein",
}

#: Whole-word alternation of the misspellings above, longest-first so that
#: overlapping keys prefer the most specific correction.  Case-insensitive so it
#: also corrects the (mixed-case) name kept for unmapped products.
_TYPO_RE = re.compile(
    r"\b(?:%s)\b"
    % "|".join(re.escape(t) for t in sorted(_TYPO_CORRECTIONS, key=len, reverse=True)),
    re.IGNORECASE,
)


def _apply_typos(text: str) -> str:
    """Correct known misspellings (whole-word, case-insensitive)."""
    return _TYPO_RE.sub(lambda m: _TYPO_CORRECTIONS[m.group(0).lower()], text)


def _normalize_product(name: str) -> str:
    """Normalise a raw protein product into a match key for a protein map."""
    # Treat separators (incl. those joining compound names) as spaces so that
    # tokens like "rdrp" in "CP/RdRp fusion" or "...; RdRp" are matchable.
    text = re.sub(r"[-_/\\();,]", " ", name.lower())
    text = re.sub(r"\s+", " ", text).strip()
    text = _QUALIFIER_RE.sub("", text)
    text = _TRAILING_RE.sub("", text)
    text = text.strip()
    # Fix known misspellings so the corrected form matches the protein map.
    text = _apply_typos(text)
    # Singularise "polymerases" for matching only, so "...RNA polymerases" reaches
    # the "...polymerase" rules while the emitted (unmapped) name keeps its plural.
    text = re.sub(r"\bpolymerases\b", "polymerase", text)
    return text


def _clean_and_capitalize(name: str) -> str:
    """Strip special characters and capitalise a lower-case leading letter."""
    for char in _SPECIAL_CHARS:
        name = name.replace(char, "")
    name = re.sub(r"\s+", " ", name).strip()
    if name and name[0].islower():
        name = name[0].upper() + name[1:]
    return name


#: A protein-map lookup: given the normalised match key and the record's
#: ``Molecule_type``, return a canonical name or ``None`` when nothing matches.
ProteinMapper = Callable[[str, str], Optional[str]]


def _standardize(name: str, mol_type: str, mapper: Optional[ProteinMapper]) -> str:
    """Run the shared cleaning pipeline, optionally applying a target ``mapper``.

    The name is cleaned (leading ``CDS:``/``ORF:`` directive stripped,
    molecular-weight and misspelling normalisation) and, when ``mapper`` is
    given, looked up against a target-specific protein map.  A map hit yields
    its canonical name; otherwise the cleaned name is kept.  Either way the
    result has special characters removed and its leading letter capitalised, and
    a name that is only a directive (``CDS``/``ORF``) or empty becomes
    ``"Unknown"``.

    Parameters
    ----------
    name : str
        Raw protein product (e.g. a FASTA-header product).
    mol_type : str
        The record's ``Molecule_type``, passed through to ``mapper`` for scoping.
    mapper : ProteinMapper, optional
        Target-specific canonicalisation; ``None`` for generic cleaning only.

    Returns
    -------
    str
        The standardised, cleaned and capitalised protein name.
    """
    # Remove any leaked NCBI "[key=value]" metadata tag (e.g. "[organism=...]").
    stripped = strip_bracket_tags(name)
    # Drop a leading CDS:/ORF: directive, so "CDS: capsid protein" still matches
    # and a bare "CDS:"/"ORF:" collapses to nothing.
    stripped = _LEADING_DIRECTIVE_RE.sub("", stripped)
    stripped = normalize_molecular_weight(stripped)
    # Correct misspellings on the kept name too, so typo variants of an unmapped
    # protein still converge (e.g. "membran protein" -> "Membrane protein").
    stripped = _apply_typos(stripped)
    # Drop leading hedging qualifiers ("putative ", "predicted ", ...) from the
    # emitted name, not just the match key.
    stripped = _strip_qualifiers(stripped)
    # Any "hypothetical ..." product (all misspellings normalised above) is
    # uninformative -> flag as Unknown so it is dropped from the FASTA and CSV.
    if stripped.strip().lower().startswith("hypothetical"):
        return "Unknown"
    # Reduce an "NSxx protein/peptide" designation to the bare "NSxx".
    stripped = _strip_designation_suffix(stripped)

    if mapper is not None:
        suggested = mapper(_normalize_product(stripped), mol_type)
        if suggested is not None:
            return _clean_and_capitalize(suggested)

    # Unmapped: a name that is only a directive (or empty) is unknown.
    result = _clean_and_capitalize(stripped)
    if result.lower() in _UNKNOWN_TOKENS:
        return "Unknown"
    return result


# ---------------------------------------------------------------------------
# Virus target: bundled viral protein map
# ---------------------------------------------------------------------------

#: Bundled viral protein-name standardization map (see the file header for the
#: normalisation / scoping rules it encodes).
_PROTEIN_MAP_TABLE = Path(__file__).resolve().parent / "data" / "viral_proteins.tsv"


def _load_protein_map() -> tuple[dict, list]:
    """Load the viral protein-name map as ``(exact, contains)`` lookups."""
    exact: dict[str, list] = {}
    contains: list = []
    if _PROTEIN_MAP_TABLE.is_file():
        with open(_PROTEIN_MAP_TABLE) as file:
            for line in file:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4 or parts[0] == "suggested_name":
                    continue
                suggested, current, match_type, scope = parts[:4]
                if match_type == "exact":
                    exact.setdefault(current, []).append((suggested, scope))
                elif match_type == "contains":
                    contains.append((current, suggested, scope))
    return exact, contains


_PROTEIN_EXACT, _PROTEIN_CONTAINS = _load_protein_map()


def _in_scope(mol_type: str, scope: str) -> bool:
    """Whether ``mol_type`` satisfies a ``molecule_type_scope`` expression."""
    mol = mol_type or ""
    for token in (part.strip() for part in scope.split(";")):
        if token == "any":
            return True
        if token == "RT" and "RT" in mol:
            return True
        if token == "RNA" and "RNA" in mol and "RT" not in mol:
            return True
        if token == "+ssRNA" and mol.startswith("ssRNA(+)"):
            return True
        if token == "-ssRNA" and mol.startswith("ssRNA(-)"):
            return True
        if token == "dsRNA" and mol.startswith("dsRNA"):
            return True
        if token == "dsDNA" and mol.startswith("dsDNA"):
            return True
        if token == "ssDNA" and mol.startswith("ssDNA"):
            return True
    return False


def _viral_mapper(normalized: str, mol_type: str) -> Optional[str]:
    """Look ``normalized`` up in the viral map, honouring the molecule-type scope."""
    for candidate, scope in _PROTEIN_EXACT.get(normalized, []):
        if _in_scope(mol_type, scope):
            return candidate
    for current, candidate, scope in _PROTEIN_CONTAINS:
        if re.search(rf"\b{re.escape(current)}\b", normalized) and _in_scope(
            mol_type, scope
        ):
            return candidate
    return None


# ---------------------------------------------------------------------------
# Per-target standardisers + dispatcher
# ---------------------------------------------------------------------------


def _standardize_virus(name: str, mol_type: str = "") -> str:
    """Standardise a viral protein name via the bundled viral protein map."""
    return _standardize(name, mol_type, _viral_mapper)


def _standardize_bacteria(name: str, mol_type: str = "") -> str:
    """Standardise a bacterial protein name (generic cleaning only, for now).

    There is no bacterial protein map yet; this is the extension point for
    bacteria-specific canonicalisation.
    """
    return _standardize(name, mol_type, None)


def _standardize_host(name: str, mol_type: str = "") -> str:
    """Standardise a host protein name (generic cleaning only).

    Host baits are gene/protein names kept as-is aside from generic cleaning;
    this is the extension point for host-specific rules.
    """
    return _standardize(name, mol_type, None)


#: Registry of per-target standardisers, keyed by ``get-databases`` target type.
_STANDARDIZERS: dict[str, Callable[[str, str], str]] = {
    "virus": _standardize_virus,
    "bacteria": _standardize_bacteria,
    "host": _standardize_host,
}


def standardize_protein(name: str, mol_type: str = "", target: str = "virus") -> str:
    """Standardise a raw protein name for the given ``target`` database.

    Parameters
    ----------
    name : str
        Raw protein product (e.g. a FASTA-header product).
    mol_type : str
        The record's ``Molecule_type``, used to scope the map rules.
    target : str
        The database target — one of ``"virus"``, ``"bacteria"`` or ``"host"``;
        selects the target-specific standardisation logic.

    Returns
    -------
    str
        The standardised, cleaned and capitalised protein name.

    Raises
    ------
    ValueError
        If ``target`` is not a known database target.
    """
    try:
        standardizer = _STANDARDIZERS[target]
    except KeyError:
        raise ValueError(f"Unknown target type: {target!r}")
    return standardizer(name, mol_type)
