"""Shared pytest fixtures for the EEfinder test suite.

Two flavours of test live in ``tests/``:

* Unit tests exercise the pure data-processing classes/functions with small,
  synthetic inputs written to ``tmp_path`` (no external binaries required).
* Integration tests (``test_integration.py``) run the ``eefinder`` CLI
  end-to-end against the files in ``test_files/`` and are skipped when the
  required binaries (blastx/diamond/bedtools) are not on ``PATH``.
"""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_FILES = REPO_ROOT / "test_files"


def pytest_addoption(parser):
    """Register ``--update-test`` to refresh the golden expected results."""
    parser.addoption(
        "--update-test",
        action="store_true",
        default=False,
        help=(
            "Overwrite test_files/expected_results/ with the outputs produced "
            "by the integration run instead of comparing against them. Use when "
            "an intended change (e.g. a dependency version bump) shifts the "
            "outputs, then commit the refreshed golden files."
        ),
    )


@pytest.fixture(scope="session")
def update_expected(request) -> bool:
    """Whether ``--update-test`` was passed (regenerate golden outputs)."""
    return bool(request.config.getoption("--update-test"))


# --------------------------------------------------------------------------- #
# Paths to the committed example inputs (used by the integration tests).
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="session")
def genome_file() -> Path:
    return TEST_FILES / "Ae_aeg_Aag2_ctg_1913.fasta"


@pytest.fixture(scope="session")
def virus_db() -> Path:
    return TEST_FILES / "virus_subset.fa"


@pytest.fixture(scope="session")
def virus_metadata() -> Path:
    return TEST_FILES / "virus_subset.csv"


@pytest.fixture(scope="session")
def filter_db() -> Path:
    return TEST_FILES / "filter_subset.fa"


@pytest.fixture(scope="session")
def expected_results() -> Path:
    """Base directory of golden main outputs, one subdir per translation method.

    ``expected_results/default/`` holds the six-frame ``blastx`` run and
    ``expected_results/gv-rv/`` the ``--translation_method gv-rv`` run; tests
    pick the subdir for the method they exercise.
    """
    return TEST_FILES / "expected_results"


# --------------------------------------------------------------------------- #
# Small synthetic helpers for the unit tests.
# --------------------------------------------------------------------------- #
def write_fasta(path: Path, records: dict[str, str]) -> Path:
    """Write ``{header: sequence}`` pairs to a FASTA file and return the path."""
    with open(path, "w") as handle:
        for header, sequence in records.items():
            handle.write(f">{header}\n{sequence}\n")
    return path


@pytest.fixture
def fasta_factory(tmp_path):
    """Return a callable that writes a FASTA file into ``tmp_path``."""

    def _make(name: str, records: dict[str, str]) -> Path:
        return write_fasta(tmp_path / name, records)

    return _make


@pytest.fixture
def blast_outfmt6(tmp_path) -> Path:
    """A minimal BLAST outfmt-6 tabular result (no header row).

    Columns: qseqid sseqid pident length mismatch gapopen qstart qend
             sstart send evalue bitscore

    Contents are crafted to exercise the filtering logic in
    :class:`eefinder.filter_table.FilterTable`:

    * ``ctg1`` has two overlapping hits within the same 100 nt window -> the
      lower-bitscore one is dropped as redundant.
    * ``ctg2`` is a negative-sense hit (qstart > qend) -> coordinates swapped.
    * ``ctg3`` is shorter than 33 aa -> removed by the length filter.
    """
    rows = [
        # qseqid sseqid       pident len mm go qstart qend sstart send evalue  bits
        [
            "ctg1",
            "PROT_A",
            "80.0",
            "100",
            "1",
            "0",
            "10",
            "310",
            "1",
            "100",
            "1e-50",
            "200.0",
        ],
        [
            "ctg1",
            "PROT_B",
            "70.0",
            "90",
            "2",
            "0",
            "40",
            "310",
            "1",
            "90",
            "1e-40",
            "150.0",
        ],
        [
            "ctg2",
            "PROT_C",
            "60.0",
            "80",
            "3",
            "0",
            "500",
            "260",
            "1",
            "80",
            "1e-30",
            "120.0",
        ],
        [
            "ctg3",
            "PROT_D",
            "90.0",
            "20",
            "0",
            "0",
            "5",
            "65",
            "1",
            "20",
            "1e-10",
            "50.0",
        ],
    ]
    path = tmp_path / "hits.blastx"
    with open(path, "w") as handle:
        for row in rows:
            handle.write("\t".join(row) + "\n")
    return path


@pytest.fixture
def taxonomy_csv(tmp_path) -> Path:
    """Metadata CSV mirroring the schema of ``test_files/virus_subset.csv``."""
    content = textwrap.dedent("""\
        Accession,Species,Genus,Family,Molecule_type,Protein,Host
        PROT_A,Species alpha,Genusalpha,Familyalpha,ssRNA(+),polyprotein,Aedes
        PROT_C,Species gamma,Genusgamma,Familygamma,ssRNA(-),glycoprotein,Culex
        """)
    path = tmp_path / "metadata.csv"
    path.write_text(content)
    return path


def binaries_available(*names: str) -> bool:
    """True only if every named executable is resolvable on ``PATH``."""
    return all(shutil.which(name) is not None for name in names)
