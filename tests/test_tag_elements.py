"""Unit tests for eefinder.tag_elements."""

from __future__ import annotations

import pandas as pd

from eefinder.tag_elements import TagElements, _list_to_string


def test_list_to_string():
    assert _list_to_string(["a", "b"]) == "a,b"
    assert _list_to_string([]) == ""


def _write_tax(path):
    pd.DataFrame(
        {
            "Element-ID": ["ctg1:100-200", "ctg1:250-300", "ctg2:100-200"],
            "Sense": ["pos", "pos", "pos"],
            "Protein-IDs": ["P1|30.0", "P2|40.0 | P3|50.0", "P4|20.0"],
            "Protein-Products": ["prot", "prot", "prot"],
            "Molecule_type": ["ssRNA", "ssRNA", "ssRNA"],
            "Family": ["FamA", "FamB", "FamC"],
            "Genus": ["GenA", "GenB", "GenC"],
            "Species": ["spA", "spB", "spC"],
            "Host": ["h", "h", "h"],
        }
    ).to_csv(path, sep="\t", index=False)


def test_tag_elements_flags_overlaps(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(tax)

    TagElements(str(tax))

    df = pd.read_csv(tax, sep="\t").set_index("Element-ID")
    # The two ctg1 elements sit within 100 nt of each other with different
    # families -> both flagged "overlaped"; the lone ctg2 element is "unique".
    assert df.loc["ctg1:100-200", "tag"] == "overlaped"
    assert df.loc["ctg1:250-300", "tag"] == "overlaped"
    assert df.loc["ctg2:100-200", "tag"] == "unique"
    assert "ctg1:250-300" in df.loc["ctg1:100-200", "Overlaped_Element_ID"]


def test_tag_elements_average_pident(tmp_path):
    tax = tmp_path / "eves.tax"
    _write_tax(tax)

    TagElements(str(tax))

    df = pd.read_csv(tax, sep="\t").set_index("Element-ID")
    assert df.loc["ctg1:100-200", "Average_pident"] == 30.0
    assert df.loc["ctg1:250-300", "Average_pident"] == 45.0  # mean(40.0, 50.0)
    assert df.loc["ctg2:100-200", "Average_pident"] == 20.0
