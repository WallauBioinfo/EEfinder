"""Unit tests for eefinder.gff.WriteGFF3."""

from __future__ import annotations

import pandas as pd

from eefinder.gff import WriteGFF3


def _write_tax(path, rows):
    columns = [
        "Element-ID",
        "Sense",
        "Protein-IDs",
        "Protein-Products",
        "Molecule_type",
        "Family",
        "Genus",
        "Species",
        "Host",
        "Overlaped_Element_ID",
        "tag",
        "Average_pident",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _read_features(path):
    lines = path.read_text().splitlines()
    header, *features = lines
    return header, [line.split("\t") for line in features]


def test_write_gff3_columns_and_coordinates(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:100-200",
                "pos",
                "P1|80.0",
                "polyprotein",
                "ssRNA(+)",
                "FamA",
                "GenA",
                "SpA",
                "Aedes",
                "",
                "unique",
                80.0,
            ],
            [
                "ctg-x:5-40",
                "neg",
                "P2|50.0 | P3|60.0",
                "glyco",
                "ssRNA(-)",
                "FamB",
                "GenB",
                "SpB",
                "Culex",
                "ctg1:100-200",
                "overlaped",
                55.0,
            ],
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    header, features = _read_features(out)
    assert header == "##gff-version 3"
    assert len(features) == 2
    by_seqid = {cols[0]: cols for cols in features}

    # seqid, source, type, start (0-based+1), end, score, strand, phase.
    # Default analysis is viral -> endogenous_viral_element.
    assert by_seqid["ctg1"][:8] == [
        "ctg1",
        "EEfinder",
        "endogenous_viral_element",
        "101",
        "200",
        "80.0",
        "+",
        ".",
    ]

    # rpartition keeps the "-" in the contig name; negative sense -> "-".
    assert by_seqid["ctg-x"][3:5] == ["6", "40"]
    assert by_seqid["ctg-x"][6] == "-"


def test_write_gff3_features_sorted_by_start(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:200-300",
                "pos",
                "P|1.0",
                "p",
                "m",
                "F",
                "G",
                "S",
                "H",
                "",
                "u",
                1.0,
            ],
            ["ctg1:50-100", "pos", "P|1.0", "p", "m", "F", "G", "S", "H", "", "u", 1.0],
            [
                "ctg1:100-150",
                "pos",
                "P|1.0",
                "p",
                "m",
                "F",
                "G",
                "S",
                "H",
                "",
                "u",
                1.0,
            ],
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    starts = [int(cols[3]) for cols in features]
    assert starts == sorted(starts) == [51, 101, 201]


def test_write_gff3_id_carries_prefix_to_match_fasta_headers(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [["ctg1:100-200", "pos", "P|1.0", "p", "m", "F", "G", "S", "H", "", "u", 1.0]],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), prefix="Aaeg")

    _, features = _read_features(out)
    attrs = dict(pair.split("=", 1) for pair in features[0][8].split(";"))
    # ID must equal the EEs.fa header "{prefix}/{Element-ID}".
    assert attrs["ID"] == "Aaeg/ctg1:100-200"


def test_write_gff3_attributes_and_escaping(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0",
                "polyprotein, partial",
                "ssRNA(+)",
                "FamA",
                "GenA",
                "SpA",
                "Aedes",
                "",
                "unique",
                90.0,
            ]
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    attrs = dict(pair.split("=", 1) for pair in features[0][8].split(";"))
    assert attrs["ID"] == "ctg1:0-10"
    assert attrs["Name"] == "SpA"
    assert attrs["family"] == "FamA"
    assert attrs["overlap_status"] == "unique"
    # The comma in the product must be percent-encoded per the GFF3 spec.
    assert attrs["product"] == "polyprotein%2C partial"


def test_write_gff3_missing_score_column(tmp_path):
    # A taxonomy table without Average_pident yields a "." score and no crash.
    tax = tmp_path / "eves.tax"
    pd.DataFrame(
        {
            "Element-ID": ["ctg1:10-20"],
            "Sense": ["pos"],
            "Family": ["FamA"],
        }
    ).to_csv(tax, sep="\t", index=False)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out))

    _, features = _read_features(out)
    assert features[0][5] == "."  # score column
    assert features[0][6] == "+"


def _single_element_tax(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0",
                "prot",
                "ssRNA",
                "FamA",
                "GenA",
                "SpA",
                "Aedes",
                "",
                "unique",
                90.0,
            ]
        ],
    )
    return tax


def test_write_gff3_analysis_virus_is_endogenous_viral_element(tmp_path):
    tax = _single_element_tax(tmp_path)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), analysis="virus")

    _, features = _read_features(out)
    assert features[0][2] == "endogenous_viral_element"


def test_write_gff3_analysis_bacteria_is_endogenous_bacterial_element(tmp_path):
    tax = _single_element_tax(tmp_path)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), analysis="bacteria")

    _, features = _read_features(out)
    assert features[0][2] == "endogenous_bacterial_element"


def test_write_gff3_explicit_feature_type_overrides_analysis(tmp_path):
    tax = _single_element_tax(tmp_path)
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), analysis="bacteria", feature_type="match")

    _, features = _read_features(out)
    assert features[0][2] == "match"


def test_write_gff3_custom_source_and_type(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(
        tax,
        [
            [
                "ctg1:0-10",
                "pos",
                "P1|90.0",
                "prot",
                "ssRNA",
                "FamA",
                "GenA",
                "SpA",
                "Aedes",
                "",
                "unique",
                90.0,
            ]
        ],
    )
    out = tmp_path / "eves.gff3"

    WriteGFF3(str(tax), str(out), source="MyTool", feature_type="match")

    _, features = _read_features(out)
    assert features[0][1] == "MyTool"
    assert features[0][2] == "match"
