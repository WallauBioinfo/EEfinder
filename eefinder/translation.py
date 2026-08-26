"""Protein prediction + coordinate traceback for the ``screening`` search.

The default similarity search translates the nucleotide query in the six reading
frames with ``blastx``/``diamond blastx``. This module implements the
alternative ``--translation-method`` modes, which predict proteins up front and
then align them with ``blastp``/``diamond blastp``:

* ``gv``    -- predict with `pyrodigal-gv <https://github.com/althonos/pyrodigal-gv>`_.
* ``rv``    -- predict with `pyrodigal-rv <https://github.com/LanderDC/pyrodigal-rv>`_.
* ``gv-rv`` -- run both predictors and drop redundancy with ``cd-hit`` (100%
  identity / 100% coverage).

Prediction writes, alongside the predicted-protein FASTA, a **coordinates TSV**
(``protein_id, contig, start, end, strand, tool``). After the protein-vs-protein
search, :func:`traceback` uses that TSV to map each hit's amino-acid coordinates
back to nucleotide coordinates on the original contig, emitting a table with the
exact same schema ``blastx`` would have produced — so every downstream step is
unchanged regardless of the chosen method.
"""

from __future__ import annotations

import re
import shlex
import subprocess

import pandas as pd
from Bio import SeqIO

from eefinder.filter_table import OUTFMT6_COLUMNS
from eefinder.log import logger

#: Translation methods selectable via ``--translation-method``.
TRANSLATION_METHODS = ("default", "gv", "rv", "gv-rv")

#: Prediction tools that translate on their own (everything but ``default``).
_PREDICTION_TOOLS = {
    "gv": ("gv",),
    "rv": ("rv",),
    "gv-rv": ("gv", "rv"),
}

#: Columns of the coordinates TSV emitted next to the predicted-protein FASTA.
COORDS_COLUMNS = ["protein_id", "contig", "start", "end", "strand", "tool"]


def _gene_finder(tool: str):
    """Return a metagenomic viral gene finder for ``gv`` or ``rv``."""
    if tool == "gv":
        import pyrodigal_gv

        return pyrodigal_gv.ViralGeneFinder(meta=True)
    if tool == "rv":
        import pyrodigal_rv

        return pyrodigal_rv.ViralGeneFinder(meta=True)
    raise ValueError(f"Unknown prediction tool: {tool!r}")


def predict_proteins(nt_fasta: str, tool: str, faa_out: str, coords_out: str) -> int:
    """Predict proteins from ``nt_fasta`` with ``tool`` (``gv``/``rv``).

    Writes the predicted proteins to ``faa_out`` (headers
    ``{contig}__{tool}__{index}``) and their nucleotide coordinates to
    ``coords_out`` (:data:`COORDS_COLUMNS`).

    Parameters
    ----------
    nt_fasta : str
        Nucleotide FASTA to predict genes on (contigs or extracted EEs).
    tool : str
        ``"gv"`` or ``"rv"``.
    faa_out : str
        Destination predicted-protein FASTA.
    coords_out : str
        Destination coordinates TSV.

    Returns
    -------
    int
        Number of predicted proteins written.
    """
    finder = _gene_finder(tool)
    n = 0
    with open(faa_out, "w") as faa, open(coords_out, "w") as coords:
        coords.write("\t".join(COORDS_COLUMNS) + "\n")
        for record in SeqIO.parse(nt_fasta, "fasta"):
            genes = finder.find_genes(str(record.seq))
            for index, gene in enumerate(genes, start=1):
                protein_id = f"{record.id}__{tool}__{index}"
                # Drop the trailing stop ("*"), which BLAST/DIAMOND reject.
                faa.write(f">{protein_id}\n{gene.translate(include_stop=False)}\n")
                strand = "+" if gene.strand >= 0 else "-"
                coords.write(
                    f"{protein_id}\t{record.id}\t{gene.begin}\t{gene.end}\t"
                    f"{strand}\t{tool}\n"
                )
                n += 1
    logger.debug(f"pyrodigal-{tool}: predicted {n} protein(s) from {nt_fasta}")
    return n


def cluster_proteins(faa_in: str, faa_out: str, threads: int) -> None:
    """Deduplicate ``faa_in`` with ``cd-hit`` at 100% identity / 100% coverage.

    Used by the ``gv-rv`` method to drop proteins predicted identically by both
    tools. Identical sequences have identical coordinates, so the coordinates
    TSV keeps every ``gv``+``rv`` entry and the cluster representative id (which
    ``cd-hit`` preserves) still resolves during :func:`traceback`.

    Parameters
    ----------
    ``cd-hit`` also writes a ``{faa_out}.clstr`` file describing the composition
    of every cluster; callers that need to know which sequences were absorbed
    into which representative read it with :func:`parse_cdhit_clusters`.

    Parameters
    ----------
    faa_in, faa_out : str
        Input and output protein FASTA.
    threads : int
        Threads for ``cd-hit`` (``-T``).
    """
    command = (
        f"cd-hit -i {faa_in} -o {faa_out} "
        f"-c 1.0 -aL 1.0 -aS 1.0 -d 0 -M 0 -T {int(threads)}"
    )
    logger.debug(f"cd-hit command: {command}")
    # capture rather than discard: on failure check=True raises with cd-hit's
    # own message attached.
    subprocess.run(
        shlex.split(command),
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        check=True,
    )


def parse_cdhit_clusters(clstr_path: str) -> "dict[str, list[str]]":
    """Parse a ``cd-hit`` ``.clstr`` file into ``cluster id -> [members]``.

    The representative -- the member ``cd-hit`` marks with ``*`` -- is always the
    first element of the list, so a caller can tell which sequence stands for the
    cluster and which ones it absorbed.

    Parameters
    ----------
    clstr_path : str
        Path to the ``.clstr`` file written next to a ``cd-hit`` output FASTA.

    Returns
    -------
    dict[str, list[str]]
        Keyed by the cluster number as ``cd-hit`` wrote it.
    """
    clusters: "dict[str, list[str]]" = {}
    cluster_id = ""
    members: "list[str]" = []
    representative = ""

    def _flush() -> None:
        if cluster_id and members:
            ordered = (
                [representative] + [m for m in members if m != representative]
                if representative
                else members
            )
            clusters[cluster_id] = ordered

    with open(clstr_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                _flush()
                cluster_id = line.split()[-1]
                members, representative = [], ""
                continue
            match = re.search(r">(.+?)\.\.\.", line)
            if not match:
                continue
            member = match.group(1)
            members.append(member)
            if line.endswith("*"):
                representative = member
    _flush()
    logger.debug(f"cd-hit: parsed {len(clusters)} cluster(s) from {clstr_path}")
    return clusters


def _aa_to_genomic(
    begin: int, end: int, strand: str, aa_start: int, aa_end: int
) -> "tuple[int, int]":
    """Map an amino-acid hit span to nucleotide coordinates on the contig.

    ``begin``/``end`` are the 1-based inclusive genomic bounds of the predicting
    CDS (``begin < end`` regardless of strand); ``aa_start``/``aa_end`` are the
    1-based inclusive hit positions in the protein. Returns ``(qstart, qend)``
    following the ``blastx`` convention: ``qstart < qend`` on the plus strand and
    ``qstart > qend`` on the minus strand, so the downstream filter infers the
    strand exactly as it does for a native ``blastx`` result.
    """
    if strand == "+":
        g_start = max(begin, begin + (aa_start - 1) * 3)
        g_end = min(end, begin + aa_end * 3 - 1)
        return g_start, g_end
    g_start = max(begin, end - aa_end * 3 + 1)
    g_end = min(end, end - (aa_start - 1) * 3)
    return g_end, g_start


def traceback(blastp_result: str, coords_tsv: str, out_blastx: str) -> None:
    """Rewrite a protein-vs-protein table with contig nucleotide coordinates.

    Reads the ``blastp``/``diamond blastp`` ``outfmt 6`` table (whose ``qseqid``
    is a predicted-protein id and whose ``qstart``/``qend`` are amino-acid
    positions), looks each protein up in ``coords_tsv`` and writes ``out_blastx``
    with the same 12 :data:`~eefinder.filter_table.OUTFMT6_COLUMNS` columns, but
    with ``qseqid`` set to the source contig and ``qstart``/``qend`` mapped to
    nucleotide coordinates.

    Hits whose protein id is absent from ``coords_tsv`` are dropped (with a
    warning), so a partially-missing map never corrupts downstream coordinates.
    """
    coords = pd.read_csv(coords_tsv, sep="\t").drop_duplicates(subset=["protein_id"])
    coord_map = {
        row.protein_id: (row.contig, int(row.start), int(row.end), row.strand)
        for row in coords.itertuples(index=False)
    }

    hits = pd.read_csv(blastp_result, sep="\t", header=None, names=OUTFMT6_COLUMNS)
    rows = []
    missing = 0
    for hit in hits.itertuples(index=False):
        info = coord_map.get(hit.qseqid)
        if info is None:
            missing += 1
            continue
        contig, begin, end, strand = info
        qstart, qend = _aa_to_genomic(
            begin, end, strand, int(hit.qstart), int(hit.qend)
        )
        rows.append(
            {
                **hit._asdict(),
                "qseqid": contig,
                "qstart": qstart,
                "qend": qend,
            }
        )
    if missing:
        logger.warning(
            f"{missing} hit(s) had no coordinate entry in {coords_tsv} and were "
            "dropped during traceback"
        )
    result = pd.DataFrame(rows, columns=OUTFMT6_COLUMNS)
    result.to_csv(out_blastx, sep="\t", header=False, index=False)


def predict_and_cluster(
    query_fasta: str, method: str, threads: int
) -> "tuple[str, str]":
    """Predict (and, for ``gv-rv``, cluster) proteins for ``query_fasta``.

    Parameters
    ----------
    query_fasta : str
        Nucleotide FASTA to predict on.
    method : str
        ``"gv"``, ``"rv"`` or ``"gv-rv"``.
    threads : int
        Threads for ``cd-hit`` (``gv-rv`` only).

    Returns
    -------
    tuple[str, str]
        ``(protein_fasta, coords_tsv)`` — the (clustered, for ``gv-rv``)
        predicted-protein FASTA and the combined coordinates TSV.
    """
    tools = _PREDICTION_TOOLS[method]
    faa_parts = []
    coords_parts = []
    for tool in tools:
        faa = f"{query_fasta}.{tool}.faa"
        coords = f"{query_fasta}.{tool}.coords.tsv"
        predict_proteins(query_fasta, tool, faa, coords)
        faa_parts.append(faa)
        coords_parts.append(coords)

    combined_coords = f"{query_fasta}.pred.coords.tsv"
    frames = [pd.read_csv(part, sep="\t") for part in coords_parts]
    pd.concat(frames, ignore_index=True).to_csv(combined_coords, sep="\t", index=False)

    combined_faa = f"{query_fasta}.pred.faa"
    with open(combined_faa, "w") as out:
        for part in faa_parts:
            with open(part) as handle:
                out.write(handle.read())

    if len(tools) > 1:
        clustered = f"{query_fasta}.pred.nr.faa"
        cluster_proteins(combined_faa, clustered, threads)
        logger.debug(f"cd-hit wrote non-redundant proteins to {clustered}")
        return clustered, combined_coords
    return combined_faa, combined_coords
