"""Write a GFF3 annotation of endogenous elements (EVEs).

Follows the GFF3 specification:
https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md
"""

from __future__ import annotations

import pandas as pd

#: Mandatory first line of a GFF3 file.
GFF3_VERSION = "##gff-version 3"

#: Column-3 feature type per analysis type. Endogenous elements are integrated
#: foreign sequences, described with the parallel ``endogenous_viral_element``
#: (EVE) / ``endogenous_bacterial_element`` terms (descriptive, not SO ids —
#: GFF3 permits free text in column 3).
FEATURE_TYPES = {
    "virus": "endogenous_viral_element",
    "bacteria": "endogenous_bacterial_element",
}

#: Analysis type assumed when none is given.
DEFAULT_ANALYSIS = "virus"

#: Map the pipeline's ``Sense`` values to GFF3 strand symbols.
_STRAND = {"pos": "+", "neg": "-"}

#: Characters that must be percent-encoded inside GFF3 column-9 values.
_ATTRIBUTE_ESCAPES = {
    "%": "%25",
    ";": "%3B",
    "=": "%3D",
    "&": "%26",
    ",": "%2C",
    "\t": "%09",
    "\n": "%0A",
    "\r": "%0D",
}

#: (GFF3 attribute tag, taxonomy-table column) pairs written to column 9, in
#: order. Reserved tags (``ID``/``Name``) are capitalised per the spec; the rest
#: are lower-case custom tags.
_ATTRIBUTE_COLUMNS = [
    ("Name", "Species"),
    ("family", "Family"),
    ("genus", "Genus"),
    ("species", "Species"),
    ("molecule_type", "Molecule_type"),
    ("product", "Protein-Products"),
    ("protein_ids", "Protein-IDs"),
    ("host", "Host"),
    ("overlap_status", "tag"),
]


def _escape(value: object) -> str:
    """Percent-encode the GFF3-reserved characters in an attribute value."""
    text = str(value)
    for char, code in _ATTRIBUTE_ESCAPES.items():
        text = text.replace(char, code)
    return text


class WriteGFF3:
    """Convert an EEfinder taxonomy table into a GFF3 annotation file.

    Each row of the taxonomy TSV becomes one feature, and features are sorted by
    sequence id then start. The ``Element-ID`` (``contig:start-end``) provides
    the sequence id and coordinates; the coordinates are BED-style (0-based,
    half-open, as emitted by ``bedtools merge``) and are converted to GFF3's
    1-based inclusive convention (``start + 1 .. end``). ``Average_pident`` is
    used as the feature score. Runs on instantiation.

    The ``ID`` attribute is prefixed with the run ``prefix`` so it matches the
    ``EEs.fa`` FASTA headers (``{prefix}/{Element-ID}``) for cross-referencing.

    Parameters
    ----------
    tax_file : str
        Path to the taxonomy TSV (e.g. ``{prefix}.EEs.tax.tsv``).
    output_file : str
        Path of the GFF3 file to write.
    prefix : str, optional
        Run prefix prepended to each element's ``ID`` so it matches the FASTA
        headers. Empty by default (``ID`` is then the bare ``Element-ID``).
    source : str, optional
        Value for GFF3 column 2 (the annotation source). Defaults to
        ``"EEfinder"``.
    analysis : str, optional
        ``"virus"`` or ``"bacteria"``; selects the column-3 feature type from
        :data:`FEATURE_TYPES` (``endogenous_viral_element`` vs
        ``endogenous_bacterial_element``). Defaults to :data:`DEFAULT_ANALYSIS`.
        Ignored when ``feature_type`` is given explicitly.
    feature_type : str, optional
        Explicit override for the GFF3 column-3 term. When ``None`` (the
        default) the term is derived from ``analysis``.
    """

    def __init__(
        self,
        tax_file: str,
        output_file: str,
        prefix: str = "",
        source: str = "EEfinder",
        analysis: str = DEFAULT_ANALYSIS,
        feature_type: "str | None" = None,
    ) -> None:
        if feature_type is None:
            if analysis not in FEATURE_TYPES:
                raise ValueError(f"Unknown analysis type: {analysis!r}")
            feature_type = FEATURE_TYPES[analysis]

        self.tax_file = tax_file
        self.output_file = output_file
        self.prefix = prefix
        self.source = source
        self.feature_type = feature_type

        self.write_gff3()

    @staticmethod
    def _parse_element_id(element_id: str) -> tuple[str, int, int]:
        """Split ``contig:start-end`` into ``(contig, start, end)``.

        ``rpartition`` on ``":"`` keeps contig names that contain ``-``.
        """
        contig, _, coords = element_id.rpartition(":")
        start, _, end = coords.partition("-")
        return contig, int(start), int(end)

    def _qualified_id(self, element_id: str) -> str:
        """Prefix the element id so it matches the ``EEs.fa`` FASTA header."""
        return f"{self.prefix}/{element_id}" if self.prefix else element_id

    def _attributes(self, row: pd.Series) -> str:
        """Build the GFF3 column-9 attribute string for one element."""
        pairs = [("ID", self._qualified_id(row["Element-ID"]))]
        for tag, column in _ATTRIBUTE_COLUMNS:
            if column not in row:
                continue
            value = row[column]
            if pd.isna(value) or str(value) == "":
                continue
            pairs.append((tag, value))
        return ";".join(f"{tag}={_escape(value)}" for tag, value in pairs)

    def write_gff3(self) -> None:
        """Read the taxonomy table and write features sorted by seqid/start."""
        df = pd.read_csv(self.tax_file, sep="\t")
        has_score = "Average_pident" in df.columns

        features = []
        for _, row in df.iterrows():
            contig, start, end = self._parse_element_id(row["Element-ID"])
            features.append((contig, start, end, row))
        features.sort(key=lambda feature: (feature[0], feature[1], feature[2]))

        with open(self.output_file, "w") as gff_out:
            gff_out.write(f"{GFF3_VERSION}\n")
            for contig, start, end, row in features:
                score = row["Average_pident"] if has_score else float("nan")
                columns = [
                    contig,
                    self.source,
                    self.feature_type,
                    str(start + 1),  # BED 0-based -> GFF3 1-based
                    str(end),
                    "." if pd.isna(score) else str(score),
                    _STRAND.get(row["Sense"], "."),
                    ".",  # phase: only meaningful for CDS features
                    self._attributes(row),
                ]
                gff_out.write("\t".join(columns) + "\n")
