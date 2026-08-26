"""Unit tests for the pure-python helpers in eefinder.bed.

The bedtools-backed classes (GetFasta, MergeBed, BedFlank) are covered by the
end-to-end integration test; here we only exercise the string/DataFrame logic.
"""

from __future__ import annotations

import pandas as pd

from eefinder.bed import GetAnnotBed, GetBed, RemoveAnnotation


def test_get_bed_parses_coordinates_from_headers(tmp_path):
    fasta = tmp_path / "eves.fa"
    fasta.write_text(">ctg1:10-310\nACGT\n>ctg2:5-99\nTTTT\n")

    GetBed(str(fasta))

    bed = pd.read_csv(
        f"{fasta}.bed", sep="\t", header=None, names=["chrom", "start", "end"]
    )
    assert bed.iloc[0].tolist() == ["ctg1", 10, 310]
    assert bed.iloc[1].tolist() == ["ctg2", 5, 99]


def test_remove_annotation_strips_pipe_suffix(tmp_path):
    annotated = tmp_path / "merged.bed"
    annotated.write_text("ctg1|FamA|GenA|pos\t10\t310\tPROT_A\n")

    RemoveAnnotation(str(annotated))

    out = pd.read_csv(f"{annotated}.fmt", sep="\t", header=None)
    assert out.iloc[0, 0] == "ctg1"


def test_get_annot_bed_genus_level_naming(tmp_path):
    tax = tmp_path / "hits.tax"
    pd.DataFrame(
        {
            "qseqid": ["ctg1:10-310"],
            "sseqid": ["PROT_A"],
            "pident": [80.0],
            "sense": ["pos"],
            "qstart": [10],
            "qend": [310],
            "Family": ["FamA"],
            "Genus": ["GenA"],
        }
    ).to_csv(tax, sep=",", index=False)

    GetAnnotBed(str(tax), merge_level="genus")

    bed = pd.read_csv(f"{tax}.bed", sep="\t", header=None)
    # formated_name collapses classified hits to contig|Family|Genus|sense.
    assert bed.iloc[0, 0] == "ctg1|FamA|GenA|pos"
    # The annotation column carries sseqid|sense|pident for later collapse.
    assert bed.iloc[0, 3] == "PROT_A|pos|80.0"
