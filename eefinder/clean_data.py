"""Sequence-level cleaning: length filtering and soft-mask filtering."""

from __future__ import annotations

from Bio import SeqIO

#: Bases counted as "masked" when deciding whether an EE sits in a repetitive
#: region: soft-masked (lowercase) bases plus ambiguous ``N``/``n``.
_MASKED_BASES = ("a", "t", "c", "g", "n", "N")


class RemoveShortSequences:
    """Drop sequences shorter than ``cutoff`` from a FASTA file.

    Writes ``{input_file}.fmt`` keeping only records whose length is greater
    than or equal to the cutoff. Runs on instantiation.

    Parameters
    ----------
    input_file : str
        Path to the FASTA file to filter.
    cutoff : int
        Minimum sequence length (nt) to keep, parsed from ``--length``.
    """

    def __init__(self, input_file: str, cutoff: int) -> None:
        self.input_file = input_file
        self.cutoff = int(cutoff)

        self.cut_seq()

    def cut_seq(self) -> None:
        """Stream records and write those at least ``cutoff`` bases long."""
        kept = (
            record
            for record in SeqIO.parse(self.input_file, "fasta")
            if len(record.seq) >= self.cutoff
        )
        with open(f"{self.input_file}.fmt", "w") as output_handle:
            SeqIO.write(kept, output_handle, "fasta")


class MaskClean:
    """Remove EEs that lie in soft-masked (repetitive) regions.

    A sequence is discarded when the fraction of masked bases
    (see :data:`_MASKED_BASES`) exceeds ``m_per`` percent. Writes
    ``{input_file}.cl``. Runs on instantiation.

    Parameters
    ----------
    input_file : str
        FASTA file with putative EEs.
    m_per : int
        Masked-percentage threshold, parsed from ``--mask_per``. Sequences with
        a masked fraction less than or equal to this value are kept.
    """

    def __init__(self, input_file: str, m_per: int) -> None:
        self.input_file = input_file
        self.m_per = int(m_per)

        self.mask_clean()

    @staticmethod
    def _masked_percentage(sequence: str) -> float:
        """Return the percentage of masked bases in ``sequence``."""
        if not sequence:
            return 0.0
        masked = sum(sequence.count(base) for base in _MASKED_BASES)
        return masked / len(sequence) * 100

    def mask_clean(self) -> None:
        """Write records whose masked percentage is within the threshold."""
        sequences: dict[str, str] = {}
        for seq_record in SeqIO.parse(self.input_file, "fasta"):
            sequence = str(seq_record.seq)
            sequence_id = str(seq_record.id)
            if self._masked_percentage(sequence) <= self.m_per:
                sequences.setdefault(sequence_id, sequence)

        with open(f"{self.input_file}.cl", "w") as output_file:
            for sequence_id, sequence in sequences.items():
                output_file.write(f">{sequence_id}\n{sequence}\n")
