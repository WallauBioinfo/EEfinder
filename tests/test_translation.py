"""Tests for the --translation-method prediction/traceback pipeline.

The pure-logic tests (coordinate mapping, traceback table rewrite, dispatch) run
anywhere. The tests that exercise real prediction/alignment are skipped unless
pyrodigal-gv/pyrodigal-rv and the cd-hit/makeblastdb/blastp binaries are present.
"""

from __future__ import annotations

import importlib.util
import subprocess

import pandas as pd
import pytest

from conftest import binaries_available
from eefinder import translation
from eefinder.filter_table import OUTFMT6_COLUMNS
from eefinder import similarity_analysis


def _module_available(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


PYRODIGAL = _module_available("pyrodigal_gv") and _module_available("pyrodigal_rv")

#: A diverse peptide (blastp SEG will not mask it) and its 1-codon-per-aa ORF.
_AA2CODON = {
    "A": "GCT",
    "R": "CGT",
    "N": "AAT",
    "D": "GAT",
    "C": "TGT",
    "Q": "CAA",
    "E": "GAA",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "L": "CTT",
    "K": "AAA",
    "M": "ATG",
    "F": "TTT",
    "P": "CCT",
    "S": "TCT",
    "T": "ACT",
    "W": "TGG",
    "Y": "TAT",
    "V": "GTT",
}
PEPTIDE = (
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"
)
ORF = "".join(_AA2CODON[a] for a in PEPTIDE) + "TAA"


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# --------------------------------------------------------------------------
# Pure coordinate logic
# --------------------------------------------------------------------------


def test_aa_to_genomic_plus_strand():
    # aa 1..10 of a CDS beginning at 100 -> nt 100..129 (ascending).
    assert translation._aa_to_genomic(100, 300, "+", 1, 10) == (100, 129)
    assert translation._aa_to_genomic(100, 300, "+", 5, 8) == (112, 123)
    # Overrun is clamped to the CDS end.
    assert translation._aa_to_genomic(100, 129, "+", 1, 99) == (100, 129)


def test_aa_to_genomic_minus_strand():
    # Minus strand counts down from `end`; emitted qstart > qend (blastx neg).
    qstart, qend = translation._aa_to_genomic(500, 800, "-", 1, 10)
    assert qstart == 800 and qend == 771
    assert qstart > qend
    # A hit deeper into the protein moves toward `begin`.
    assert translation._aa_to_genomic(500, 800, "-", 5, 8) == (788, 777)


def test_traceback_rewrites_query_coordinates(tmp_path):
    coords = tmp_path / "coords.tsv"
    coords.write_text(
        "protein_id\tcontig\tstart\tend\tstrand\ttool\n"
        "p1\tc1\t100\t300\t+\tgv\n"
        "p2\tc1\t500\t800\t-\tgv\n"
    )
    # outfmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    blastp = tmp_path / "hits.blastp"
    blastp.write_text(
        "p1\trefA\t95.0\t10\t0\t0\t1\t10\t1\t10\t1e-20\t80.0\n"
        "p2\trefB\t90.0\t10\t0\t0\t1\t10\t1\t10\t1e-18\t70.0\n"
        "pX\trefC\t99.0\t10\t0\t0\t1\t10\t1\t10\t1e-30\t99.0\n"  # missing -> dropped
    )
    out = tmp_path / "out.blastx"
    translation.traceback(str(blastp), str(coords), str(out))

    result = pd.read_csv(out, sep="\t", header=None, names=OUTFMT6_COLUMNS)
    assert list(result["qseqid"]) == ["c1", "c1"]  # pX dropped, contig substituted
    p1 = result.iloc[0]
    assert (int(p1.qstart), int(p1.qend)) == (100, 129)  # plus
    p2 = result.iloc[1]
    assert (int(p2.qstart), int(p2.qend)) == (800, 771)  # minus (qstart > qend)
    # Non-coordinate columns are carried through unchanged.
    assert float(p1.pident) == 95.0 and float(p2.bitscore) == 70.0


# --------------------------------------------------------------------------
# Dispatch: one translation_method controls the backend
# --------------------------------------------------------------------------


def test_similarity_search_default_uses_blastx(monkeypatch):
    calls = {}
    monkeypatch.setattr(
        similarity_analysis,
        "runblastx",
        lambda q, d, t: calls.setdefault("blastx", (q, d, t)),
    )
    monkeypatch.setattr(
        similarity_analysis,
        "run_predicted_search",
        lambda *a: calls.setdefault("predicted", a),
    )
    similarity_analysis.SimilaritySearch("q.fa", "db.fa", 2, "blastx", "default")
    assert "blastx" in calls and "predicted" not in calls


@pytest.mark.parametrize("method", ["gv", "rv", "gv-rv"])
def test_similarity_search_prediction_methods_route_to_predicted(monkeypatch, method):
    calls = {}
    monkeypatch.setattr(
        similarity_analysis,
        "runblastx",
        lambda *a: calls.setdefault("blastx", a),
    )
    monkeypatch.setattr(
        similarity_analysis,
        "run_predicted_search",
        lambda q, d, t, m, meth: calls.setdefault("predicted", meth),
    )
    similarity_analysis.SimilaritySearch("q.fa", "db.fa", 2, "blastx", method)
    assert calls.get("predicted") == method and "blastx" not in calls


# --------------------------------------------------------------------------
# Real prediction / alignment (skipped without the tools)
# --------------------------------------------------------------------------


@pytest.mark.skipif(not PYRODIGAL, reason="requires pyrodigal-gv and pyrodigal-rv")
@pytest.mark.parametrize("tool", ["gv", "rv"])
def test_predict_proteins_writes_faa_and_coords(tmp_path, tool):
    genome = tmp_path / "genome.fa"
    genome.write_text(f">c1\n{'N' * 60 + ORF + 'N' * 30}\n")
    faa = tmp_path / "p.faa"
    coords = tmp_path / "p.tsv"
    n = translation.predict_proteins(str(genome), tool, str(faa), str(coords))
    assert n >= 1
    assert "*" not in faa.read_text()  # trailing stop stripped
    df = pd.read_csv(coords, sep="\t")
    assert list(df.columns) == translation.COORDS_COLUMNS
    row = df.iloc[0]
    assert row.contig == "c1" and row.tool == tool and row.strand in ("+", "-")
    assert 1 <= int(row.start) < int(row.end) <= len(genome.read_text())


@pytest.mark.skipif(
    not (PYRODIGAL and binaries_available("cd-hit")),
    reason="requires pyrodigal + cd-hit",
)
def test_predict_and_cluster_gv_rv_dedups(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(f">c1\n{'N' * 60 + ORF + 'N' * 30}\n")
    faa, coords = translation.predict_and_cluster(str(genome), "gv-rv", 1)

    combined = tmp_path / "genome.fa.pred.faa"
    n_before = sum(1 for line in combined.open() if line.startswith(">"))
    n_after = sum(1 for line in open(faa) if line.startswith(">"))
    assert n_before == 2  # gv + rv predicted the same protein
    assert n_after == 1  # cd-hit collapsed the identical pair
    # The combined coords TSV keeps both tools' entries with identical coords.
    df = pd.read_csv(coords, sep="\t")
    assert set(df["tool"]) == {"gv", "rv"}
    assert df["start"].nunique() == 1 and df["end"].nunique() == 1


@pytest.mark.skipif(
    not (PYRODIGAL and binaries_available("makeblastdb", "blastp")),
    reason="requires pyrodigal + makeblastdb/blastp",
)
@pytest.mark.parametrize("strand", ["plus", "minus"])
def test_run_predicted_search_maps_back_to_nucleotides(tmp_path, strand):
    contig_seq = "N" * 60 + (ORF if strand == "plus" else _revcomp(ORF)) + "N" * 30
    genome = tmp_path / "genome.fa"
    genome.write_text(f">contig1\n{contig_seq}\n")

    # Protein DB = the reference peptide itself, so blastp finds the predicted ORF.
    db = tmp_path / "db.faa"
    db.write_text(f">ref\n{PEPTIDE}\n")
    subprocess.run(
        ["makeblastdb", "-dbtype", "prot", "-in", str(db)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )

    similarity_analysis.run_predicted_search(str(genome), str(db), 1, "blastx", "gv")
    out = genome.parent / "genome.fa.blastx"
    hits = pd.read_csv(out, sep="\t", header=None, names=OUTFMT6_COLUMNS)
    assert len(hits) >= 1
    hit = hits.iloc[0]
    assert hit.qseqid == "contig1"  # protein id mapped back to the contig
    lo, hi = sorted((int(hit.qstart), int(hit.qend)))
    assert 61 <= lo and hi <= len(contig_seq)  # within the ORF region on the contig
    if strand == "plus":
        assert int(hit.qstart) < int(hit.qend)
    else:
        assert int(hit.qstart) > int(hit.qend)  # blastx neg convention
