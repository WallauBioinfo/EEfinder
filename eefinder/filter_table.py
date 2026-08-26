"""Filter a translated-search tabular result down to non-redundant hits."""

from __future__ import annotations

import glob
import os
import shutil

import pandas as pd

#: Standard BLAST/DIAMOND ``outfmt 6`` columns, in order.
OUTFMT6_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

#: Columns of the filtered table, i.e. outfmt6 plus EEfinder annotations.
FILTERED_COLUMNS = OUTFMT6_COLUMNS + ["sense", "bed_name", "tag"]

#: Minimum alignment length (aa) for a hit to be retained.
MIN_HIT_LENGTH = 33

#: Rows processed per chunk when streaming very large result tables.
CHUNK_SIZE = 200_000


class FilterTable:
    """Collapse redundant translated-search hits into a filtered table.

    For each query, hits falling in the same ``rangejunction``-sized window on
    the same strand are deduplicated (the highest-bitscore hit wins, since the
    table is pre-sorted by bitscore). Negative-strand hits have their
    coordinates normalised so ``qstart < qend``. Hits shorter than
    :data:`MIN_HIT_LENGTH` are dropped.

    Writes ``{blast_result}.filtred`` and, for EE searches, the
    ``{blast_result}.filtred.bed`` coordinate file. Runs on instantiation.

    Parameters
    ----------
    blast_result : str
        Path to the raw ``outfmt 6`` table.
    rangejunction : int
        Window size (nt) used to merge redundant hits, from ``--range_junction``.
    tag : str
        ``"EE"`` for the endogenous-element search or ``"HOST"`` for the
        host-gene bait search; controls how ``bed_name`` is built.
    out_dir : str
        Output directory, used for scratch chunk files.
    """

    def __init__(
        self, blast_result: str, rangejunction: int, tag: str, out_dir: str
    ) -> None:
        self.blast_result = blast_result
        self.rangejunction = rangejunction
        self.tag = tag
        self.out_dir = out_dir

        self.filter_blast()

    def _annotate_chunk(self, df: pd.DataFrame) -> pd.DataFrame:
        """Assign strand, normalise coordinates and build ``bed_name``/``tag``."""
        df["sense"] = df["sense"].astype(object)
        df.loc[df["qstart"].astype(int) > df["qend"].astype(int), "sense"] = "neg"
        df.loc[df["qend"].astype(int) > df["qstart"].astype(int), "sense"] = "pos"

        # Normalise negative-strand hits so start < end.
        neg = df["sense"] == "neg"
        df.loc[neg, ["qstart", "qend"]] = df.loc[neg, ["qend", "qstart"]].values

        if self.tag == "EE":
            df["tag"] = "EE"
            df["bed_name"] = df.apply(
                lambda x: f"{x['qseqid']}:{x['qstart']}-{x['qend']}", axis=1
            )
        else:
            df["tag"] = "HOST"
            df["bed_name"] = df["qseqid"]

        df["evalue"] = pd.to_numeric(df["evalue"], downcast="float")
        df = df[df.length >= MIN_HIT_LENGTH]
        return df[FILTERED_COLUMNS]

    def filter_blast(self) -> None:
        """Run the full filter: annotate, deduplicate and write outputs."""
        df = pd.read_csv(
            self.blast_result, sep="\t", header=None, names=OUTFMT6_COLUMNS
        ).sort_values(by="bitscore", ascending=False)
        df["sense"] = ""
        df["bed_name"] = ""
        df["tag"] = ""
        df.to_csv(f"{self.blast_result}.csv", sep="\t")

        tmp_path = f"{self.out_dir}/tmp/"
        os.makedirs(tmp_path, exist_ok=True)

        chunks = pd.read_csv(f"{self.blast_result}.csv", sep="\t", chunksize=CHUNK_SIZE)
        for count, chunk in enumerate(chunks):
            annotated = self._annotate_chunk(chunk)
            annotated.to_csv(f"{tmp_path}chunk.{count}.tsv", sep="\t", index=False)

        filtered = pd.concat(
            (pd.read_csv(chunk, sep="\t") for chunk in glob.glob(f"{tmp_path}/*.tsv")),
            ignore_index=True,
        )
        # Merge hits that share a strand and fall in the same coordinate window.
        filtered["qstart_rng"] = filtered.qstart.floordiv(self.rangejunction)
        filtered["qend_rng"] = filtered.qend.floordiv(self.rangejunction)
        filtered = filtered.drop_duplicates(
            subset=["qseqid", "qstart_rng", "sense"]
        ).sort_values(by=["qseqid"])

        filtered.to_csv(
            f"{self.blast_result}.filtred",
            sep="\t",
            index=False,
            columns=FILTERED_COLUMNS,
        )
        if self.tag == "EE":
            filtered.to_csv(
                f"{self.blast_result}.filtred.bed",
                header=False,
                sep="\t",
                index=False,
                columns=["qseqid", "qstart", "qend"],
            )
        shutil.rmtree(tmp_path, ignore_errors=True)
