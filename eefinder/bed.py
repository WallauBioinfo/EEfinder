"""BED construction and bedtools wrappers used across the pipeline.

The bedtools-backed steps (:class:`GetFasta`, :class:`MergeBed`,
:class:`BedFlank`) shell out to the ``bedtools`` binary; the remaining classes
manipulate coordinates with pandas / string parsing only.
"""

from __future__ import annotations

import re
import shlex
import subprocess

import numpy as np
import pandas as pd


class GetFasta:
    """Extract sequences for BED intervals via ``bedtools getfasta``.

    Runs on instantiation.

    Parameters
    ----------
    input_file : str
        FASTA to extract from (parsed from ``--genome_file`` upstream).
    bed_file : str
        BED file of intervals to extract.
    out_file : str
        Output FASTA path.
    """

    def __init__(self, input_file: str, bed_file: str, out_file: str) -> None:
        self.input_file = input_file
        self.bed_file = bed_file
        self.out_file = out_file

        self.get_fasta()

    def get_fasta(self) -> None:
        """Run ``bedtools getfasta``."""
        cmd = (
            f"bedtools getfasta "
            f"-fi {self.input_file} "
            f"-bed {self.bed_file} "
            f"-fo {self.out_file}"
        )
        subprocess.run(
            shlex.split(cmd),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


class GetAnnotBed:
    """Build an annotated BED so truncated EEs of a taxon can be merged.

    Elements are named ``contig|Family|Genus|sense`` (genus level) or
    ``contig|Family|sense`` (family level) so that ``bedtools merge`` joins only
    neighbouring hits of the same taxon and strand. Writes
    ``{blast_tax_info}.bed``. Runs on instantiation.

    Parameters
    ----------
    blast_tax_info : str
        Taxonomy signature CSV from :class:`~eefinder.get_taxonomy.GetTaxonomy`.
    merge_level : str
        ``"genus"`` or ``"family"``; the level at which to merge, from
        ``--merge_level``.
    """

    def __init__(self, blast_tax_info: str, merge_level: str) -> None:
        self.blast_tax_info = blast_tax_info
        self.merge_level = merge_level

        self.get_annotated_bed()

    def get_annotated_bed(self) -> None:
        """Compute the merge-friendly element names and write the BED file."""
        df = pd.read_csv(self.blast_tax_info, sep=",")
        df["qseqid"] = df["qseqid"].str.replace(r"\:.*", "", regex=True)
        df["sseqid"] = df["sseqid"] + "|" + df["sense"] + "|" + df["pident"].astype(str)
        df["Family"] = df["Family"].fillna("Unknown")
        df["Genus"] = df["Genus"].fillna("Unknown")

        if self.merge_level == "genus":
            df["formated_name"] = np.where(
                df["Genus"] != "Unknown",
                df["qseqid"]
                + "|"
                + df["Family"]
                + "|"
                + df["Genus"]
                + "|"
                + df["sense"],
                df["qseqid"] + "|" + df["sseqid"] + "|" + df["Genus"],
            )
        else:
            df["formated_name"] = np.where(
                df["Family"] != "Unknown",
                df["qseqid"] + "|" + df["Family"] + "|" + df["sense"],
                df["qseqid"] + "|" + df["sseqid"] + "|" + df["Family"],
            )

        bed = df[["formated_name", "qstart", "qend", "sseqid"]].copy()
        bed = bed.sort_values(["formated_name", "qstart"], ascending=(True, True))
        bed.to_csv(f"{self.blast_tax_info}.bed", index=False, header=False, sep="\t")


class RemoveAnnotation:
    """Strip the ``|Family|Genus|sense`` annotation from merged BED names.

    Writes ``{bed_annotated_merged_file}.fmt``. Runs on instantiation.

    Parameters
    ----------
    bed_annotated_merged_file : str
        Merged BED file produced by :class:`MergeBed`.
    """

    def __init__(self, bed_annotated_merged_file: str) -> None:
        self.bed_annotated_merged_file = bed_annotated_merged_file

        self.reformat_bed()

    def reformat_bed(self) -> None:
        """Reduce the first column back to the bare contig name."""
        df = pd.read_csv(self.bed_annotated_merged_file, sep="\t", header=None)
        df.iloc[:, 0] = df.iloc[:, 0].str.replace(r"\|.*", "", regex=True)
        df.to_csv(
            f"{self.bed_annotated_merged_file}.fmt",
            index=False,
            header=False,
            sep="\t",
        )


class MergeBed:
    """Merge nearby intervals of the same annotated name via ``bedtools merge``.

    Writes ``{bed_annotated_file}.merge``. Runs on instantiation.

    Parameters
    ----------
    bed_annotated_file : str
        Annotated BED from :class:`GetAnnotBed`.
    limit_merge : int
        Maximum gap (nt) between intervals to merge, from ``--limit``.
    """

    def __init__(self, bed_annotated_file: str, limit_merge: int) -> None:
        self.bed_annotated_file = bed_annotated_file
        self.limit_merge = limit_merge

        self.merge_bed()

    def merge_bed(self) -> None:
        """Run ``bedtools merge`` collapsing the annotation column."""
        cmd = (
            f"bedtools merge "
            f"-d {int(self.limit_merge)} "
            f"-i {self.bed_annotated_file} "
            f"-c 4 "
            f"-o collapse "
            f'-delim " AND "'
        )
        with open(f"{self.bed_annotated_file}.merge", "w") as merge_output:
            subprocess.run(shlex.split(cmd), stdout=merge_output)


class BedFlank:
    """Extend intervals by a flank on each side via ``bedtools slop``.

    Writes ``{input_file}.flank``. Runs on instantiation.

    Parameters
    ----------
    input_file : str
        BED file from :class:`GetBed`.
    length_file : str
        Genome length index from :class:`~eefinder.get_length.GetLength`.
    flank_region : int
        Number of bases to add on each side, from ``--flank``.
    """

    def __init__(self, input_file: str, length_file: str, flank_region: int) -> None:
        self.input_file = input_file
        self.length_file = length_file
        self.flank_region = flank_region

        self.bedtools_flank()

    def bedtools_flank(self) -> None:
        """Run ``bedtools slop`` to add the flanking regions."""
        cmd = (
            f"bedtools slop "
            f"-i {self.input_file} "
            f"-g {self.length_file} "
            f"-b {int(self.flank_region)}"
        )
        with open(f"{self.input_file}.flank", "w") as flank_out:
            subprocess.run(shlex.split(cmd), stdout=flank_out)


class GetBed:
    """Derive a BED file from ``>contig:start-end`` FASTA headers.

    Writes ``{input_file}.bed``. Runs on instantiation.

    Parameters
    ----------
    input_file : str
        FASTA whose headers encode ``contig:start-end`` coordinates.
    """

    def __init__(self, input_file: str) -> None:
        self.input_file = input_file

        self.get_bed()

    def get_bed(self) -> None:
        """Parse each header into ``contig<TAB>start<TAB>end``."""
        with open(self.input_file) as fasta_in:
            headers = [
                line[1:].rstrip("\n") for line in fasta_in if line.startswith(">")
            ]

        with open(f"{self.input_file}.bed", "w") as bed_out:
            for header in headers:
                # rpartition on ":" so contig names containing "-" are kept whole.
                contig, _, coords = header.rpartition(":")
                start, _, end = coords.partition("-")
                bed_out.write(f"{contig}\t{start}\t{end}\n")
