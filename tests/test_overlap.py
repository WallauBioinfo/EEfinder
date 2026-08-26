"""Unit tests for eefinder.overlap."""

from __future__ import annotations

import pandas as pd

from eefinder.overlap import FilterOverlap, _element_length, elements_to_remove

TAX_COLUMNS = [
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


def _tax_frame(rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=TAX_COLUMNS)


def _row(element_id, family, tag, partners=""):
    return [
        element_id,
        "pos",
        "P|1.0",
        "prot",
        "ssRNA",
        family,
        "Gen",
        "Sp",
        "Host",
        partners,
        tag,
        1.0,
    ]


def test_element_length_parses_coordinates():
    assert _element_length("ctg1:100-250") == 150
    # Contig names with extra separators keep only the final coordinate span.
    assert _element_length("ctg-x:5-40") == 35


def test_keep_removes_nothing():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
        ]
    )
    assert elements_to_remove(df, "keep", []) == set()


def test_targets_drops_non_target_members_of_a_cluster_with_a_target():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
            _row("ctg2:0-500", "FamC", "unique"),
        ]
    )
    # The cluster contains a FamA (target) member, so its FamB member is
    # dropped; the unique FamC element is untouched (never overlaped).
    assert elements_to_remove(df, "targets", ["FamA"]) == {"ctg1:150-400"}


def test_targets_keeps_whole_cluster_without_a_target_family():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamB", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamC", "overlaped", "ctg1:100-200"),
        ]
    )
    # No member of the cluster is a target family, so the targets logic does not
    # apply and the entire cluster is kept.
    assert elements_to_remove(df, "targets", ["FamA"]) == set()


def test_targets_keeps_cluster_where_all_members_are_targets():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
        ]
    )
    # Both families are targets, so nothing is dropped.
    assert elements_to_remove(df, "targets", ["FamA", "FamB"]) == set()


def test_targets_applies_per_cluster_independently():
    df = _tax_frame(
        [
            # Cluster 1 (ctg1): has a FamA target -> drop the FamB member.
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
            # Cluster 2 (ctg2): no target member -> keep both.
            _row("ctg2:100-200", "FamB", "overlaped", "ctg2:150-400"),
            _row("ctg2:150-400", "FamC", "overlaped", "ctg2:100-200"),
        ]
    )
    assert elements_to_remove(df, "targets", ["FamA"]) == {"ctg1:150-400"}


def test_targets_chained_cluster_drops_all_non_targets():
    df = _tax_frame(
        [
            # A-B-C chain (A-C do not directly overlap) is a single cluster.
            _row("ctg1:0-1000", "FamA", "overlaped", "ctg1:100-200"),
            _row("ctg1:100-200", "FamB", "overlaped", "ctg1:0-1000,ctg1:150-400"),
            _row("ctg1:150-400", "FamC", "overlaped", "ctg1:100-200"),
        ]
    )
    # FamA is in the cluster, so both non-target members (FamB, FamC) are dropped.
    assert elements_to_remove(df, "targets", ["FamA"]) == {
        "ctg1:100-200",
        "ctg1:150-400",
    }


def test_non_targets_drops_listed_families_from_cluster():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
        ]
    )
    # FamB is a non-target (drop-list): drop it, keep the rest.
    assert elements_to_remove(df, "targets", [], ["FamB"]) == {"ctg1:150-400"}


def test_non_targets_never_wipes_a_fully_listed_cluster():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
        ]
    )
    # Every member is on the drop-list, so the cluster is kept untouched.
    assert elements_to_remove(df, "targets", [], ["FamA", "FamB"]) == set()


def test_non_targets_leaves_clusters_without_a_listed_family():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
        ]
    )
    # FamC is not present in the cluster, so nothing is dropped.
    assert elements_to_remove(df, "targets", [], ["FamC"]) == set()


def test_longest_drops_shorter_of_each_overlap():
    df = _tax_frame(
        [
            _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),  # len 100
            _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),  # len 250
        ]
    )
    assert elements_to_remove(df, "longest", []) == {"ctg1:100-200"}


def test_longest_keeps_when_element_is_the_longest_partner():
    df = _tax_frame(
        [
            _row("ctg1:0-1000", "FamA", "overlaped", "ctg1:100-200,ctg1:300-350"),
            _row("ctg1:100-200", "FamB", "overlaped", "ctg1:0-1000"),
            _row("ctg1:300-350", "FamC", "overlaped", "ctg1:0-1000"),
        ]
    )
    assert elements_to_remove(df, "longest", []) == {"ctg1:100-200", "ctg1:300-350"}


def _write_tax(path, rows):
    _tax_frame(rows).to_csv(path, sep="\t", index=False)


def _write_fasta(path, headers):
    with open(path, "w") as handle:
        for header in headers:
            handle.write(f">PFX/{header}\nACGT\n")


def test_filter_overlap_splits_kept_and_removed(tmp_path):
    tax = tmp_path / "eves.tax"
    fasta = tmp_path / "eves.fa"
    removed_tax = tmp_path / "removed.tax"
    removed_fasta = tmp_path / "removed.fa"

    rows = [
        _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
        _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
    ]
    _write_tax(tax, rows)
    _write_fasta(fasta, ["ctg1:100-200", "ctg1:150-400"])

    FilterOverlap(
        str(fasta),
        str(tax),
        "targets",
        ["FamA"],
        str(removed_fasta),
        str(removed_tax),
    )

    kept = pd.read_csv(tax, sep="\t")
    removed = pd.read_csv(removed_tax, sep="\t")
    assert list(kept["Element-ID"]) == ["ctg1:100-200"]
    assert list(removed["Element-ID"]) == ["ctg1:150-400"]

    assert ">PFX/ctg1:100-200" in fasta.read_text()
    assert ">PFX/ctg1:150-400" not in fasta.read_text()
    assert ">PFX/ctg1:150-400" in removed_fasta.read_text()


def test_filter_overlap_uses_non_target_families(tmp_path):
    tax = tmp_path / "eves.tax"
    fasta = tmp_path / "eves.fa"
    removed_tax = tmp_path / "removed.tax"
    removed_fasta = tmp_path / "removed.fa"

    rows = [
        _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
        _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
    ]
    _write_tax(tax, rows)
    _write_fasta(fasta, ["ctg1:100-200", "ctg1:150-400"])

    FilterOverlap(
        str(fasta),
        str(tax),
        "targets",
        [],
        str(removed_fasta),
        str(removed_tax),
        non_target_families=["FamB"],
    )

    kept = pd.read_csv(tax, sep="\t")
    removed = pd.read_csv(removed_tax, sep="\t")
    assert list(kept["Element-ID"]) == ["ctg1:100-200"]
    assert list(removed["Element-ID"]) == ["ctg1:150-400"]


def test_filter_overlap_keep_strategy_is_a_no_op(tmp_path):
    tax = tmp_path / "eves.tax"
    fasta = tmp_path / "eves.fa"
    rows = [
        _row("ctg1:100-200", "FamA", "overlaped", "ctg1:150-400"),
        _row("ctg1:150-400", "FamB", "overlaped", "ctg1:100-200"),
    ]
    _write_tax(tax, rows)
    _write_fasta(fasta, ["ctg1:100-200", "ctg1:150-400"])

    FilterOverlap(
        str(fasta),
        str(tax),
        "keep",
        [],
        str(tmp_path / "removed.fa"),
        str(tmp_path / "removed.tax"),
    )

    kept = pd.read_csv(tax, sep="\t")
    assert list(kept["Element-ID"]) == ["ctg1:100-200", "ctg1:150-400"]
    assert (tmp_path / "removed.tax").exists()  # empty removed table still written
