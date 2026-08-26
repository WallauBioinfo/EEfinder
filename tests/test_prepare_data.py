"""Unit tests for eefinder.prepare_data."""

from __future__ import annotations

from Bio import SeqIO

from eefinder.clean_data import RemoveShortSequences
from eefinder.prepare_data import InsertPrefix, PrepareGenome


def test_insert_prefix_rewrites_all_headers(fasta_factory, tmp_path):
    fasta = fasta_factory("genome.fasta", {"ctg1": "ACGTACGT", "ctg2": "TTTTGGGG"})

    InsertPrefix(str(fasta), "PFX", str(tmp_path))

    out = tmp_path / "PFX.rn"
    headers = [rec.id for rec in SeqIO.parse(str(out), "fasta")]
    assert headers == ["PFX/ctg1", "PFX/ctg2"]


def test_insert_prefix_preserves_sequences(fasta_factory, tmp_path):
    fasta = fasta_factory("genome.fasta", {"ctg1": "ACGTACGT"})

    InsertPrefix(str(fasta), "PFX", str(tmp_path))

    out = tmp_path / "PFX.rn"
    record = next(SeqIO.parse(str(out), "fasta"))
    assert str(record.seq) == "ACGTACGT"


class TestPrepareGenome:
    """Prefixing and length filtering happen in one pass, writing one file."""

    def _genome(self, tmp_path):
        genome = tmp_path / "genome.fasta"
        genome.write_text(
            ">ctg1 first contig\n" + "A" * 120 + "\n"
            ">ctg2 short one\n" + "C" * 30 + "\n"
            ">ctg3\n" + "G" * 200 + "\n"
        )
        return genome

    def test_only_the_prepared_file_is_written(self, tmp_path):
        """The intermediate `.rn` doubled the disk cost of a genome."""
        outdir = tmp_path / "out"
        outdir.mkdir()
        PrepareGenome(str(self._genome(tmp_path)), "Aaeg", str(outdir), 100)

        assert [p.name for p in outdir.iterdir()] == ["Aaeg.rn.fmt"]

    def test_headers_are_prefixed_and_short_contigs_dropped(self, tmp_path):
        outdir = tmp_path / "out"
        outdir.mkdir()
        step = PrepareGenome(str(self._genome(tmp_path)), "Aaeg", str(outdir), 100)

        records = list(SeqIO.parse(outdir / "Aaeg.rn.fmt", "fasta"))
        assert [record.id for record in records] == ["Aaeg/ctg1", "Aaeg/ctg3"]
        # The whole original header is kept, prefixed once.
        assert records[0].description == "Aaeg/ctg1 first contig"
        assert (step.total, step.kept) == (3, 2)

    def test_it_matches_the_two_step_path_it_replaces(self, tmp_path):
        """Byte-for-byte, so no downstream result can shift."""
        genome = self._genome(tmp_path)
        old, new = tmp_path / "old", tmp_path / "new"
        old.mkdir()
        new.mkdir()

        InsertPrefix(str(genome), "Aaeg", str(old))
        RemoveShortSequences(str(old / "Aaeg.rn"), 100)
        PrepareGenome(str(genome), "Aaeg", str(new), 100)

        assert (old / "Aaeg.rn.fmt").read_bytes() == (new / "Aaeg.rn.fmt").read_bytes()

    def test_a_cutoff_of_zero_keeps_everything(self, tmp_path):
        outdir = tmp_path / "out"
        outdir.mkdir()
        step = PrepareGenome(str(self._genome(tmp_path)), "Aaeg", str(outdir), 0)
        assert step.kept == 3
