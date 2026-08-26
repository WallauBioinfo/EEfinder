"""Unit tests for eefinder.clean_data."""

from __future__ import annotations

from Bio import SeqIO

from eefinder.clean_data import MaskClean, RemoveShortSequences


def test_remove_short_sequences_drops_below_cutoff(fasta_factory):
    fasta = fasta_factory(
        "seqs.fasta",
        {"long": "A" * 50, "short": "A" * 5},
    )

    RemoveShortSequences(str(fasta), cutoff=10)

    kept = {rec.id for rec in SeqIO.parse(f"{fasta}.fmt", "fasta")}
    assert kept == {"long"}


def test_remove_short_sequences_keeps_equal_to_cutoff(fasta_factory):
    fasta = fasta_factory("seqs.fasta", {"exact": "A" * 10})

    RemoveShortSequences(str(fasta), cutoff=10)

    kept = {rec.id for rec in SeqIO.parse(f"{fasta}.fmt", "fasta")}
    assert kept == {"exact"}


def test_mask_clean_removes_soft_masked_sequences(fasta_factory):
    # lowercase bases + N count as masked; the fully-lowercase record is 100%
    # masked and must be dropped at the default 50% threshold.
    fasta = fasta_factory(
        "eves.fasta",
        {"clean": "ACGTACGTAC", "masked": "acgtacgtac"},
    )

    MaskClean(str(fasta), m_per=50)

    kept = {rec.id for rec in SeqIO.parse(f"{fasta}.cl", "fasta")}
    assert kept == {"clean"}


def test_mask_clean_threshold_is_inclusive(fasta_factory):
    # 5/10 lowercase == exactly 50%, which is <= threshold -> kept.
    fasta = fasta_factory("eves.fasta", {"half": "ACGTAacgta"})

    MaskClean(str(fasta), m_per=50)

    kept = {rec.id for rec in SeqIO.parse(f"{fasta}.cl", "fasta")}
    assert kept == {"half"}
