"""Unit tests for taxonomy expansion (no network: the CLI is faked)."""

import json
import subprocess

import pytest

from eefinder import taxon_exclusion
from eefinder.taxon_exclusion import (
    Expansion,
    TaxonNode,
    batch_taxa,
    expand_taxon_excluding,
    resolve_taxon,
    summarize_taxa,
)

# A small tree:  1 -> 10 -> {100, 200}, 100 -> {1000, 2000}
TREE = {
    1: {"name": "Root", "parents": [], "children": [10]},
    10: {"name": "Mid", "parents": [1], "children": [100, 200]},
    100: {"name": "Left", "parents": [1, 10], "children": [1000, 2000]},
    200: {"name": "Right", "parents": [1, 10], "children": []},
    1000: {"name": "Target", "parents": [1, 10, 100], "children": []},
    2000: {"name": "Sibling", "parents": [1, 10, 100], "children": []},
    9: {"name": "Outside", "parents": [], "children": []},
}
BY_NAME = {node["name"].lower(): tax_id for tax_id, node in TREE.items()}


def _fake_datasets(monkeypatch, calls=None):
    """Replace subprocess.run with a lookup into TREE."""

    def fake_run(args, **kwargs):
        wanted = args[args.index("taxon") + 1].split(",")
        if calls is not None:
            calls.append(tuple(wanted))
        lines = []
        for item in wanted:
            tax_id = BY_NAME.get(item.lower())
            if tax_id is None:
                try:
                    tax_id = int(item)
                except ValueError:
                    tax_id = None
            if tax_id is None or tax_id not in TREE:
                return subprocess.CompletedProcess(
                    args, 1, "", f"Error: '{item}' is not recognized."
                )
            node = TREE[tax_id]
            lines.append(
                json.dumps(
                    {
                        "taxonomy": {
                            "tax_id": tax_id,
                            "current_scientific_name": {"name": node["name"]},
                            "parents": node["parents"],
                            "children": node["children"],
                        }
                    }
                )
            )
        return subprocess.CompletedProcess(args, 0, "\n".join(lines) + "\n", "")

    monkeypatch.setattr(taxon_exclusion.subprocess, "run", fake_run)


def test_summarize_taxa_parses_nodes(monkeypatch):
    _fake_datasets(monkeypatch)
    nodes = summarize_taxa(["10", "100"])
    assert set(nodes) == {10, 100}
    assert nodes[10] == TaxonNode(10, "Mid", (1,), (100, 200))


def test_summarize_taxa_empty_input_makes_no_call(monkeypatch):
    _fake_datasets(monkeypatch)
    assert summarize_taxa([]) == {}


def test_summarize_taxa_raises_on_unknown(monkeypatch):
    _fake_datasets(monkeypatch)
    with pytest.raises(RuntimeError, match="could not look up taxonomy"):
        summarize_taxa(["nope"])


def test_resolve_taxon_by_name(monkeypatch):
    _fake_datasets(monkeypatch)
    assert resolve_taxon("Target").tax_id == 1000


def test_no_exclusions_returns_the_root_untouched(monkeypatch):
    _fake_datasets(monkeypatch)
    result = expand_taxon_excluding("10", [])
    assert result.taxa == (10,)
    assert not result.pruned


def test_excluding_a_leaf_keeps_its_siblings(monkeypatch):
    _fake_datasets(monkeypatch)
    result = expand_taxon_excluding("10", ["1000"])
    # 10 -> {100, 200}; 100 -> {1000, 2000}. Dropping 1000 keeps 200 and 2000.
    assert set(result.taxa) == {200, 2000}
    assert result.pruned
    assert [node.tax_id for node in result.excluded] == [1000]


def test_excluding_a_whole_branch_drops_it(monkeypatch):
    _fake_datasets(monkeypatch)
    result = expand_taxon_excluding("10", ["100"])
    assert set(result.taxa) == {200}


def test_exclusion_by_scientific_name(monkeypatch):
    _fake_datasets(monkeypatch)
    assert set(expand_taxon_excluding("10", ["Target"]).taxa) == {200, 2000}


def test_two_exclusions_compose(monkeypatch):
    _fake_datasets(monkeypatch)
    result = expand_taxon_excluding("10", ["1000", "200"])
    assert set(result.taxa) == {2000}
    assert {node.tax_id for node in result.excluded} == {1000, 200}


def test_taxon_outside_the_root_is_reported_not_raised(monkeypatch):
    _fake_datasets(monkeypatch)
    result = expand_taxon_excluding("100", ["200"])
    assert result.taxa == (100,)
    assert not result.pruned
    assert [node.tax_id for node in result.not_below_root] == [200]


def test_excluding_the_root_itself_is_an_error(monkeypatch):
    _fake_datasets(monkeypatch)
    with pytest.raises(RuntimeError, match="nothing would be left"):
        expand_taxon_excluding("10", ["10"])


def test_unknown_root_is_an_error(monkeypatch):
    _fake_datasets(monkeypatch)
    with pytest.raises(RuntimeError):
        expand_taxon_excluding("nope", ["1000"])


def test_lookups_are_cached_per_expansion(monkeypatch):
    calls = []
    _fake_datasets(monkeypatch, calls)
    expand_taxon_excluding("10", ["1000", "2000"])
    # Each node is fetched at most once even though both exclusions walk 10->100.
    flat = [item for call in calls for item in call]
    assert len(flat) == len(set(flat))


def test_batch_taxa_splits_at_the_inputfile_limit():
    assert batch_taxa(list(range(250)), limit=100) == [
        list(range(100)),
        list(range(100, 200)),
        list(range(200, 250)),
    ]
    assert batch_taxa([1, 2], limit=100) == [[1, 2]]


def test_expansion_pruned_flag():
    root = TaxonNode(1, "Root", (), ())
    assert not Expansion((1,), root, (), ()).pruned
    assert Expansion((2,), root, (root,), ()).pruned
