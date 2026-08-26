"""Attach taxonomy metadata to filtered hits and build the EE taxonomy table."""

from __future__ import annotations

import csv
import re

import pandas as pd
from Bio import SeqIO

#: Column order of the taxonomy tables written for the user.
TAXONOMY_COLUMNS = [
    "Element-ID",
    "Sense",
    "Protein-IDs",
    "Protein-Products",
    "Molecule_type",
    "Family",
    "Genus",
    "Species",
    "Host",
]

# Column indices in the taxonomy-signature CSV produced by GetTaxonomy
# (the filtered outfmt6 columns followed by the joined metadata columns).
_ACCESSION_COL = 1
_SPECIES_COL = 15
_GENUS_COL = 16
_FAMILY_COL = 17
_MOLTYPE_COL = 18
_PRODUCT_COL = 19
_HOST_COL = 20


class GetTaxonomy:
    """Join filtered hits with the metadata CSV to build a taxonomy signature.

    Left-merges the deduplicated BLAST table on ``sseqid`` with the metadata
    (whose ``Accession`` column is treated as ``sseqid``) and writes
    ``{blast_file}.tax``. Runs on instantiation.

    Parameters
    ----------
    blast_file : str
        TSV of filtered BLAST results (must contain an ``sseqid`` column).
    tax_file : str
        Metadata CSV, parsed from ``--dbmetadata``.
    """

    def __init__(self, blast_file: str, tax_file: str) -> None:
        self.blast_file = blast_file
        self.tax_file = tax_file

        self.get_taxonomy()

    def get_taxonomy(self) -> None:
        """Left-merge hits with metadata and write ``{blast_file}.tax``."""
        df_blast_file = pd.read_csv(self.blast_file, sep="\t")
        df_tax_file = pd.read_csv(self.tax_file)
        df_tax_file.rename(columns={"Accession": "sseqid"}, inplace=True)
        df_merged = pd.merge(df_blast_file, df_tax_file, on="sseqid", how="left")
        df_merged.to_csv(f"{self.blast_file}.tax", index=False, header=True)


class GetFinalTaxonomy:
    """Assemble a per-element taxonomy row from the merged BED file.

    For each merged element, the constituent protein accessions are looked up
    in the taxonomy signature and their products/molecule types/families/etc.
    are aggregated (multiple hits are joined with ``" AND "``). Writes
    ``{bed_formated}.fa.tax``. Runs on instantiation.

    Parameters
    ----------
    bed_formated : str
        Merged BED file (contig, start, end, collapsed protein annotation).
    taxonomy_info : str
        Taxonomy signature CSV produced by :class:`GetTaxonomy`.
    """

    def __init__(self, bed_formated: str, taxonomy_info: str) -> None:
        self.bed_formated = bed_formated
        self.taxonomy_info = taxonomy_info

        self.get_final_taxonomy()

    def _load_taxonomy_rows(self) -> list[list[str]]:
        """Read the taxonomy signature CSV once into a list of rows."""
        with open(self.taxonomy_info) as prot_info:
            return list(csv.reader(prot_info, delimiter=","))

    def get_final_taxonomy(self) -> None:
        """Build one taxonomy row per merged element and write the table."""
        taxonomy_rows = self._load_taxonomy_rows()

        with open(self.bed_formated) as bed_merge_file:
            output_rows = [
                self._build_row(line, taxonomy_rows)
                for line in csv.reader(bed_merge_file, delimiter="\t")
            ]

        with open(f"{self.bed_formated}.fa.tax", "w") as bed_merge_tax_out:
            writer = csv.writer(bed_merge_tax_out, delimiter="\t")
            writer.writerow(TAXONOMY_COLUMNS)
            writer.writerows(output_rows)

    def _build_row(self, line: list[str], taxonomy_rows: list[list[str]]) -> list[str]:
        """Resolve a single merged element into a taxonomy row."""
        element_merged_id = f"{line[0]}:{line[1]}-{line[2]}"

        sense = ""
        if "pos" in line[3]:
            sense = "pos"
            line[3] = re.sub(r"\|pos", "", line[3])
        elif "neg" in line[3]:
            sense = "neg"
            line[3] = re.sub(r"\|neg", "", line[3])
        protein_ids = line[3]

        protein_terms = ""
        mol_type = ""
        family = ""
        genus = ""
        species = ""
        host = ""

        if "AND" in protein_ids:
            # Multiple proteins were merged into this element.
            protein_ids = re.sub("AND", "|", line[3])
            for prot in taxonomy_rows:
                if prot[_ACCESSION_COL] not in protein_ids:
                    continue
                if prot[_PRODUCT_COL] not in protein_terms:
                    protein_terms += prot[_PRODUCT_COL] + " AND "
                    mol_type = prot[_MOLTYPE_COL]
                    family = prot[_FAMILY_COL]
                if prot[_GENUS_COL] not in genus:
                    genus += prot[_GENUS_COL] + " AND "
                if prot[_SPECIES_COL] not in species:
                    species += prot[_SPECIES_COL] + " AND "
                if prot[_HOST_COL] not in host:
                    host += prot[_HOST_COL] + " AND "
        else:
            for prot in taxonomy_rows:
                if prot[_ACCESSION_COL] in protein_ids:
                    protein_terms = prot[_PRODUCT_COL]
                    mol_type = prot[_MOLTYPE_COL]
                    family = prot[_FAMILY_COL]
                    genus = prot[_GENUS_COL]
                    species = prot[_SPECIES_COL]
                    host = prot[_HOST_COL]

        protein_terms = re.sub(r" AND $", "", protein_terms)
        genus = re.sub(r" AND $", "", genus)
        species = re.sub(r" AND $", "", species)
        host = re.sub(r" AND $", "", host)

        family = family or "Unclassified"
        genus = genus or "Unclassified"
        species = species or "Unclassified"
        host = host or "Undefined"

        return [
            element_merged_id,
            sense,
            protein_ids,
            protein_terms,
            mol_type,
            family,
            genus,
            species,
            host,
        ]


class GetCleanedTaxonomy:
    """Subset the taxonomy table to the elements surviving mask-cleaning.

    Reads the cleaned EE FASTA and keeps only the taxonomy rows whose
    ``Element-ID`` matches a retained record. Writes ``{cleaned_file}.tax``.
    Runs on instantiation. Parsed by the ``--clean_masked`` option.

    Parameters
    ----------
    cleaned_file : str
        Cleaned EE FASTA produced by :class:`~eefinder.clean_data.MaskClean`.
    taxonomy_file : str
        The full taxonomy table (TSV) to subset.
    """

    def __init__(self, cleaned_file: str, taxonomy_file: str) -> None:
        self.cleaned_file = cleaned_file
        self.taxonomy_file = taxonomy_file

        self.get_cleaned_taxonomy()

    def get_cleaned_taxonomy(self) -> None:
        """Keep taxonomy rows whose Element-ID is present in the cleaned FASTA."""
        with open(self.taxonomy_file) as tax_file:
            rows = list(csv.reader(tax_file, delimiter="\t"))
        if not rows:
            return
        header, data_rows = rows[0], rows[1:]

        # The taxonomy table stores Element-IDs with the "PREFIX/" removed (see
        # TagElements), so strip it from the FASTA ids before matching. Keeping
        # the source header ensures the column count matches the copied rows.
        kept_ids = {
            re.sub(r".*/", "", seq_record.id)
            for seq_record in SeqIO.parse(self.cleaned_file, "fasta")
        }
        output_rows = [header]
        output_rows.extend(row for row in data_rows if row and row[0] in kept_ids)

        with open(f"{self.cleaned_file}.tax", "w") as output_file:
            csv.writer(output_file, delimiter="\t").writerows(output_rows)
