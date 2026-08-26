"""Unit tests for eefinder.get_taxonomy and eefinder.compare_results."""

from __future__ import annotations

import csv

import pandas as pd

from eefinder.compare_results import CompareResults
from eefinder.get_taxonomy import GetCleanedTaxonomy, GetTaxonomy


def test_get_taxonomy_left_merges_metadata(tmp_path, taxonomy_csv):
    blast = tmp_path / "hits.tsv"
    pd.DataFrame(
        {
            "qseqid": ["ctg1:1-100", "ctg2:1-100"],
            "sseqid": ["PROT_A", "PROT_Z"],  # PROT_Z is absent from the metadata
            "pident": [80.0, 60.0],
        }
    ).to_csv(blast, sep="\t", index=False)

    GetTaxonomy(str(blast), str(taxonomy_csv))

    merged = pd.read_csv(f"{blast}.tax")
    assert "Family" in merged.columns
    row_a = merged[merged["sseqid"] == "PROT_A"].iloc[0]
    assert row_a["Family"] == "Familyalpha"
    # Left join keeps the unmatched hit with empty taxonomy.
    assert merged[merged["sseqid"] == "PROT_Z"]["Family"].isna().all()


def test_compare_results_host_hit_removes_element(tmp_path):
    vir = tmp_path / "vir.filtred"
    host = tmp_path / "host.filtred"

    pd.DataFrame(
        {
            "qseqid": ["x", "y"],
            "sseqid": ["V1", "V2"],
            "bitscore": [200.0, 100.0],
            "bed_name": ["ee1:1-10", "ee2:1-10"],
            "tag": ["EE", "EE"],
        }
    ).to_csv(vir, sep="\t", index=False)

    pd.DataFrame(
        {
            "qseqid": ["ee1:1-10"],
            "sseqid": ["H1"],
            "bitscore": [300.0],  # beats the viral hit for ee1
            "bed_name": ["ee1:1-10"],
            "tag": ["HOST"],
        }
    ).to_csv(host, sep="\t", index=False)

    CompareResults(str(vir), str(host))

    nr = pd.read_csv(f"{host}.concat.nr", sep="\t")
    # ee1 loses to the higher-bitscore host hit and is filtered out; only the
    # viral-only ee2 survives.
    assert set(nr["qseqid"]) == {"ee2:1-10"}
    assert (nr["tag"] == "EE").all()


def test_get_cleaned_taxonomy_matches_ids_ignoring_prefix(tmp_path):
    # The main taxonomy stores Element-IDs without a "PREFIX/"; the cleaned
    # FASTA keeps it. GetCleanedTaxonomy must still match, and preserve the
    # source header (so column counts stay consistent for downstream tagging).
    tax = tmp_path / "main.tax"
    with open(tax, "w") as handle:
        handle.write("Element-ID\tFamily\ttag\n")
        handle.write("ctg:1-10\tFamA\tunique\n")
        handle.write("ctg:20-30\tFamB\tunique\n")

    cleaned = tmp_path / "clean.fa"
    cleaned.write_text(">PFX/ctg:1-10\nACGT\n")  # only the first element survives

    GetCleanedTaxonomy(str(cleaned), str(tax))

    with open(f"{cleaned}.tax") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    assert rows[0] == ["Element-ID", "Family", "tag"]  # header preserved
    assert [row[0] for row in rows[1:]] == ["ctg:1-10"]
