"""Input preparation: tag genome FASTA headers and drop the short contigs.

:class:`PrepareGenome` is what the pipeline runs: it does both in a single pass,
writing only ``{prefix}.rn.fmt``. The two operations used to be separate steps
chained through a ``{prefix}.rn`` file, which meant writing the whole genome to
disk twice -- 3.2 GB of intermediates for a 1.6 GB genome, of which the first
copy was never read again after the second step.

:class:`InsertPrefix` and :class:`RemoveShortSequences` remain available for
callers that need one operation on its own.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from eefinder.log import logger


class InsertPrefix:
    """Prefix every FASTA header of a genome so element IDs are traceable.

    Each header ``>contig`` becomes ``>{prefix}/contig``. The class runs on
    instantiation and writes ``{outdir}/{prefix}.rn``.

    Parameters
    ----------
    input_file : str
        Path to the input genome FASTA file (nucleotides).
    prefix : str
        Prefix to prepend to each sequence header.
    outdir : str
        Output directory for the prefixed file.

    Example
    -------
    >>> InsertPrefix("genome.fasta", "Aaeg", "results")  # doctest: +SKIP
    """

    def __init__(self, input_file: str, prefix: str, outdir: str) -> None:
        self.input_file = input_file
        self.prefix = prefix
        self.outdir = outdir

        self.insert_prefix()

    def insert_prefix(self) -> None:
        """Stream the FASTA, prepending ``{prefix}/`` to each header line."""
        output_file = Path(self.outdir) / f"{self.prefix}.rn"
        with open(self.input_file) as reader, open(output_file, "w") as writer:
            for line in reader:
                if line.startswith(">"):
                    line = f">{self.prefix}/{line[1:]}"
                writer.write(line)


def prefix_record(record, prefix: str):
    """Prefix a record's id and description with ``{prefix}/``.

    Both are set so that Biopython writes the header exactly as the two-step
    path did: the whole original header line, prefixed once.
    """
    record.description = f"{prefix}/{record.description}"
    record.id = f"{prefix}/{record.id}"
    return record


class PrepareGenome:
    """Prefix the headers and drop short contigs in one pass.

    Writes ``{outdir}/{prefix}.rn.fmt``, the file every later step reads. Runs on
    instantiation.

    Parameters
    ----------
    input_file : str
        Input genome FASTA (nucleotides).
    prefix : str
        Prefix prepended to every header, so element IDs stay traceable to the
        run.
    outdir : str
        Output directory.
    cutoff : int
        Minimum contig length (nt) to keep, from ``--length``.

    Attributes
    ----------
    total, kept : int
        Contigs read and contigs written.
    """

    def __init__(self, input_file: str, prefix: str, outdir: str, cutoff: int) -> None:
        self.input_file = input_file
        self.prefix = prefix
        self.outdir = outdir
        self.cutoff = int(cutoff)
        self.total = 0
        self.kept = 0

        self.prepare()

    def _records(self):
        """Yield the prefixed records that are long enough, counting both."""
        for record in SeqIO.parse(self.input_file, "fasta"):
            self.total += 1
            if len(record.seq) >= self.cutoff:
                self.kept += 1
                yield prefix_record(record, self.prefix)

    def prepare(self) -> None:
        """Stream the genome once, writing the prepared FASTA."""
        output_file = Path(self.outdir) / f"{self.prefix}.rn.fmt"
        with open(output_file, "w") as handle:
            SeqIO.write(self._records(), handle, "fasta")
        logger.debug(
            f"PrepareGenome: {self.kept} of {self.total} contig(s) kept "
            f"(>= {self.cutoff} nt) in {output_file}"
        )
