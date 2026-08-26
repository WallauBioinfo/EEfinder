"""Run the similarity search: translated (blastx) or predicted (blastp).

The ``default`` translation method runs ``blastx``/``diamond blastx`` (six-frame
translation). The ``gv``/``rv``/``gv-rv`` methods predict proteins first (see
:mod:`eefinder.translation`) and align them with ``blastp``/``diamond blastp``,
then map the amino-acid coordinates back to nucleotides so the output schema is
identical. Both similarity searches in ``screening`` (the main EE search and the
host-bait search) go through :class:`SimilaritySearch`, so a single
``translation_method`` value keeps them consistent.
"""

from __future__ import annotations

import shlex
import subprocess

from Bio.Blast.Applications import NcbiblastxCommandline

from eefinder import translation
from eefinder.log import logger

#: E-value cutoff shared by both search backends.
EVALUE_CUTOFF = 0.00001


def runblastx(query_file: str, database_file: str, threads: int) -> None:
    """Run NCBI ``blastx`` writing tabular (outfmt 6) to ``{query_file}.blastx``."""
    cline = NcbiblastxCommandline(
        query=query_file,
        db=database_file,
        out=f"{query_file}.blastx",
        outfmt=6,
        word_size=3,
        evalue=EVALUE_CUTOFF,
        num_threads=threads,
        matrix="BLOSUM45",
        max_intron_length=100,
        soft_masking="true",
    )
    cline()


def rundiamond(query_file: str, database_file: str, threads: int, mode: str) -> None:
    """Run ``diamond blastx`` in the given sensitivity ``mode``."""
    clinedmd = (
        f"diamond blastx "
        f"-p {int(threads)} "
        f"-d {database_file}.dmnd "
        f"-f 6 "
        f"-q {query_file} "
        f"-o {query_file}.blastx "
        f"-e {EVALUE_CUTOFF} "
        f"--matrix BLOSUM45 "
        f"-k 500 "
        f"--max-hsps 0 "
        f"--{mode}"
    )
    subprocess.run(
        shlex.split(clinedmd),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def runblastp(query_file: str, database_file: str, threads: int, out: str) -> None:
    """Run NCBI ``blastp`` (protein query vs protein DB) writing outfmt 6."""
    command = (
        f"blastp "
        f"-query {query_file} "
        f"-db {database_file} "
        f"-out {out} "
        f"-outfmt 6 "
        f"-evalue {EVALUE_CUTOFF} "
        f"-matrix BLOSUM45 "
        f"-num_threads {int(threads)}"
    )
    subprocess.run(shlex.split(command), check=True)


def rundiamond_blastp(
    query_file: str, database_file: str, threads: int, mode: str, out: str
) -> None:
    """Run ``diamond blastp`` in the given sensitivity ``mode`` writing outfmt 6."""
    command = (
        f"diamond blastp "
        f"-p {int(threads)} "
        f"-d {database_file}.dmnd "
        f"-f 6 "
        f"-q {query_file} "
        f"-o {out} "
        f"-e {EVALUE_CUTOFF} "
        f"--matrix BLOSUM45 "
        f"-k 500 "
        f"--max-hsps 0 "
        f"--{mode}"
    )
    subprocess.run(
        shlex.split(command),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def run_predicted_search(
    query_file: str,
    database_file: str,
    threads: int,
    mode: str,
    method: str,
) -> None:
    """Predict proteins, align with blastp, and write ``{query_file}.blastx``.

    Shared by both ``screening`` similarity searches for the ``gv``/``rv``/
    ``gv-rv`` methods: predicts proteins from the nucleotide ``query_file`` (see
    :func:`eefinder.translation.predict_and_cluster`), runs ``blastp`` or
    ``diamond blastp`` against ``database_file``, then traces the amino-acid hit
    coordinates back to contig nucleotides so the emitted ``{query_file}.blastx``
    matches a native ``blastx`` table.
    """
    protein_fasta, coords_tsv = translation.predict_and_cluster(
        query_file, method, threads
    )
    blastp_out = f"{query_file}.pred.blastp"
    logger.debug(
        f"Predicted search ({method}, {mode}): {protein_fasta} vs {database_file}"
    )
    if mode == "blastx":
        runblastp(protein_fasta, database_file, threads, blastp_out)
    else:
        rundiamond_blastp(protein_fasta, database_file, threads, mode, blastp_out)
    translation.traceback(blastp_out, coords_tsv, f"{query_file}.blastx")


class SimilaritySearch:
    """Run a BLAST or DIAMOND translated search of ``query`` vs ``database``.

    Runs on instantiation and writes tabular output to ``{query_file}.blastx``.

    Parameters
    ----------
    query_file : str
        Nucleotide FASTA used as the query.
    database_file : str
        Protein database created by :class:`~eefinder.make_database.MakeDB`.
    threads : int
        Number of threads for the search.
    mode : str
        ``"blastx"`` for NCBI BLAST, otherwise a DIAMOND sensitivity mode.
    translation_method : str
        ``"default"`` (six-frame ``blastx``/``diamond blastx``) or a prediction
        method (``"gv"``/``"rv"``/``"gv-rv"``) that predicts proteins and aligns
        them with ``blastp``/``diamond blastp`` (see :func:`run_predicted_search`).
        Both ``screening`` searches pass the same value.
    """

    def __init__(
        self,
        query_file: str,
        database_file: str,
        threads: int,
        mode: str,
        translation_method: str = "default",
    ) -> None:
        self.query_file = query_file
        self.database_file = database_file
        self.threads = threads
        self.mode = mode
        self.translation_method = translation_method

        self.similarity_search()

    def similarity_search(self) -> None:
        """Dispatch on the translation method, then the BLAST/DIAMOND backend."""
        if self.translation_method != "default":
            run_predicted_search(
                self.query_file,
                self.database_file,
                self.threads,
                self.mode,
                self.translation_method,
            )
        elif self.mode == "blastx":
            runblastx(self.query_file, self.database_file, self.threads)
        else:
            rundiamond(self.query_file, self.database_file, self.threads, self.mode)
