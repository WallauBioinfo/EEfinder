"""Unit tests for eefinder.get_length.GetLength."""

from __future__ import annotations

from eefinder.get_length import GetLength


def test_get_length_writes_id_and_length(fasta_factory):
    fasta = fasta_factory("genome.rn.fmt", {"ctg1": "A" * 40, "ctg2": "C" * 15})

    GetLength(str(fasta))

    # GetLength appends the ".rn.fmt.lenght" suffix to the input path.
    lengths = {}
    with open(f"{fasta}.rn.fmt.lenght") as handle:
        for line in handle:
            name, size = line.split()
            lengths[name] = int(size)

    assert lengths == {"ctg1": 40, "ctg2": 15}
