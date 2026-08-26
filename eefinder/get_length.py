"""Produce a genome length index (``<id>\\t<length>``) for bedtools slop."""

from __future__ import annotations

from Bio import SeqIO


class GetLength:
    """Write a two-column ``id<TAB>length`` file for every FASTA record.

    The output is used as the genome file (``-g``) for ``bedtools slop`` when
    extracting flanking regions. Runs on instantiation.

    Parameters
    ----------
    input_file : str
        Formatted genome FASTA (``{prefix}.rn.fmt``) produced upstream.

    Notes
    -----
    For backwards compatibility the output path keeps its historical
    (double-suffixed) name ``{input_file}.rn.fmt.lenght``.
    """

    def __init__(self, input_file: str) -> None:
        self.input_file = input_file

        self.get_length()

    def get_length(self) -> None:
        """Write ``id<TAB>length`` for each record in ``input_file``."""
        with open(f"{self.input_file}.rn.fmt.lenght", "w") as output_length:
            for seq_record in SeqIO.parse(self.input_file, "fasta"):
                output_length.write(f"{seq_record.id}\t{len(seq_record)}\n")
