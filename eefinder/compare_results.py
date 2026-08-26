"""Reconcile viral/bacterial EE hits against host-gene bait hits."""

from __future__ import annotations

import pandas as pd


class CompareResults:
    """Keep only putative EEs that out-score their host-bait matches.

    The two filtered BLAST tables are concatenated per query; for each query
    only the highest-bitscore hit is kept. Queries whose best hit is a host
    gene (``tag == "HOST"``) are dropped, leaving genuine EEs. Runs on
    instantiation and writes ``{host_result}.concat`` (all hits) and
    ``{host_result}.concat.nr`` (deduplicated, EE-only).

    Parameters
    ----------
    vir_result : str
        Filtered BLAST table against the EE (virus/bacteria) database.
    host_result : str
        Filtered BLAST table against the host-gene bait database.
    """

    def __init__(self, vir_result: str, host_result: str) -> None:
        self.vir_result = vir_result
        self.host_result = host_result

        self.compare_results()

    def compare_results(self) -> None:
        """Merge, deduplicate by query and retain the EE-only best hits."""
        df_vir = pd.read_csv(self.vir_result, sep="\t")
        # Align the viral queries on their coordinate-tagged bed name so they
        # collide with the host hits extracted from the same region.
        df_vir["qseqid"] = df_vir["bed_name"]
        df_host = pd.read_csv(self.host_result, sep="\t")

        df_hybrid = pd.concat([df_vir, df_host], ignore_index=True)
        df_hybrid = df_hybrid.sort_values(by=["qseqid", "bitscore"], ascending=False)
        df_hybrid.to_csv(f"{self.host_result}.concat", sep="\t", index=False)

        best_per_query = df_hybrid.drop_duplicates(subset=["qseqid"])
        ee_only = best_per_query[best_per_query.tag == "EE"]
        ee_only.to_csv(f"{self.host_result}.concat.nr", sep="\t", index=False)
