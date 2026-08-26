"""Build BLAST or DIAMOND protein databases from a FASTA file."""

from __future__ import annotations

import shlex
import subprocess

from Bio.Blast.Applications import NcbimakeblastdbCommandline

#: DIAMOND / BLAST modes that use NCBI BLAST rather than DIAMOND.
BLAST_MODES = ("blastx", "tblastn")


def makeblastdb(data: str, db_type: str) -> None:
    """Create an NCBI BLAST database from ``data`` (``nucl`` or ``prot``)."""
    clinedb = NcbimakeblastdbCommandline(dbtype=db_type, input_file=data)
    clinedb()


def makediamonddb(data: str, threads: int) -> None:
    """Create a DIAMOND database (``{data}.dmnd``) with the BLOSUM45 matrix."""
    clinedb = (
        f"diamond makedb "
        f"--db {data} "
        f"--in {data} "
        f"--threads {int(threads)} "
        f"--matrix BLOSUM45"
    )
    subprocess.run(
        shlex.split(clinedb),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


class MakeDB:
    """Create the similarity-search database for the chosen ``mode``.

    BLAST modes use ``makeblastdb``; DIAMOND modes use ``diamond makedb``.
    Runs on instantiation.

    Parameters
    ----------
    mode : str
        Search mode selected with ``--mode``.
    data : str
        Path to the database FASTA file, parsed from ``--database``.
    db_type : str
        ``"nucl"`` or ``"prot"`` (only used by BLAST).
    threads : int
        Number of threads (only used by DIAMOND).
    """

    def __init__(self, mode: str, data: str, db_type: str, threads: int) -> None:
        self.mode = mode
        self.data = data
        self.db_type = db_type
        self.threads = threads

        self.make_db()

    def make_db(self) -> None:
        """Dispatch to the BLAST or DIAMOND database builder."""
        if self.mode in BLAST_MODES:
            makeblastdb(self.data, self.db_type)
        else:
            makediamonddb(self.data, self.threads)
