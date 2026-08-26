"""Unit tests for eefinder.filter_table.FilterTable."""

from __future__ import annotations

import pandas as pd

from eefinder.filter_table import FilterTable


def test_filter_table_ee_outputs(blast_outfmt6, tmp_path):
    FilterTable(str(blast_outfmt6), rangejunction=100, tag="EE", out_dir=str(tmp_path))

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")

    # ctg3 (length 20 < 33) is dropped, and the two overlapping ctg1 hits in the
    # same 100 nt window collapse to the higher-bitscore one (PROT_A).
    assert set(filtred["qseqid"]) == {"ctg1", "ctg2"}
    assert "PROT_B" not in set(filtred["sseqid"])
    assert "PROT_D" not in set(filtred["sseqid"])


def test_filter_table_negative_strand_is_swapped(blast_outfmt6, tmp_path):
    FilterTable(str(blast_outfmt6), rangejunction=100, tag="EE", out_dir=str(tmp_path))

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")
    ctg2 = filtred[filtred["qseqid"] == "ctg2"].iloc[0]

    # Original hit was qstart=500 > qend=260 -> negative sense, coords swapped.
    assert ctg2["sense"] == "neg"
    assert ctg2["qstart"] == 260
    assert ctg2["qend"] == 500


def test_filter_table_bed_name_and_bed_file(blast_outfmt6, tmp_path):
    FilterTable(str(blast_outfmt6), rangejunction=100, tag="EE", out_dir=str(tmp_path))

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")
    ctg1 = filtred[filtred["qseqid"] == "ctg1"].iloc[0]
    assert ctg1["bed_name"] == "ctg1:10-310"
    assert (filtred["tag"] == "EE").all()

    bed = pd.read_csv(
        f"{blast_outfmt6}.filtred.bed",
        sep="\t",
        header=None,
        names=["qseqid", "qstart", "qend"],
    )
    assert len(bed) == 2
    assert list(bed.columns) == ["qseqid", "qstart", "qend"]


def test_filter_table_host_tag_uses_qseqid_as_bed_name(blast_outfmt6, tmp_path):
    FilterTable(
        str(blast_outfmt6), rangejunction=100, tag="HOST", out_dir=str(tmp_path)
    )

    filtred = pd.read_csv(f"{blast_outfmt6}.filtred", sep="\t")
    assert (filtred["tag"] == "HOST").all()
    # For HOST hits the bed_name is simply the query id, not coordinate-tagged.
    assert (filtred["bed_name"] == filtred["qseqid"]).all()
