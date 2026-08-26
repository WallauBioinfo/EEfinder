"""The download must finish when the *archive* is done, not when the tool exits.

The NCBI `datasets` client writes the package, prints its completed validation
bar and then keeps running, which froze `get-databases` on a download that had
actually succeeded. These tests stand in a fake `datasets` that behaves exactly
that way, so the symptom is reproduced without the network.
"""

from __future__ import annotations

import io
import json
import os
import time
import zipfile

import pandas as pd
import pytest

from eefinder.get_databases import GetDatabases

from conftest import binaries_available

PROTEIN = (
    ">YP_000001.1 glycoprotein [organism=Mivirus chuvi]\n"
    # Longer than cd-hit's throw-away length (-l 10), which would silently
    # discard it.
    "MKALLVGTSGAGKSTLLQALNRLYELDSGSIRIDG\n"
)
REPORT = json.dumps(
    {"virus": {"organismName": "Mivirus chuvi", "lineage": [{"name": "Chuviridae"}]}}
)


def _write_package(path: str) -> None:
    """Write a minimal but valid datasets package."""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("ncbi_dataset/data/protein.faa", PROTEIN)
        archive.writestr("ncbi_dataset/data/data_report.jsonl", REPORT + "\n")


@pytest.fixture
def hanging_datasets(tmp_path):
    """A `datasets` that delivers the package and then never exits."""
    package = tmp_path / "package.zip"
    _write_package(str(package))

    script = tmp_path / "datasets"
    script.write_text(
        "#!/bin/sh\n"
        # The real client takes --filename <path>; find it the same way.
        'while [ $# -gt 0 ]; do if [ "$1" = "--filename" ]; then out="$2"; fi;'
        " shift; done\n"
        f'cp "{package}" "$out"\n'
        'printf "Downloading: %s    104kB valid data package\\n" "$out"\n'
        'printf "Validating package files [====] 100%% 5/5\\n"\n'
        "sleep 300\n"  # <- the behaviour that used to hang the pipeline
    )
    script.chmod(0o755)
    return str(script)


@pytest.fixture
def exiting_datasets(tmp_path):
    """A well-behaved `datasets` that exits after writing the package."""
    package = tmp_path / "package.zip"
    _write_package(str(package))

    script = tmp_path / "datasets"
    script.write_text(
        "#!/bin/sh\n"
        'while [ $# -gt 0 ]; do if [ "$1" = "--filename" ]; then out="$2"; fi;'
        " shift; done\n"
        f'cp "{package}" "$out"\n'
        'printf "Validating package files [====] 100%% 5/5\\n"\n'
    )
    script.chmod(0o755)
    return str(script)


def _download(tmp_path, datasets_bin, cluster=False, **kwargs):
    return GetDatabases(
        dataset="virus",
        taxon="Chuviridae",
        outdir=str(tmp_path / "db"),
        prefix="chuvi",
        cluster=cluster,
        datasets_bin=datasets_bin,
        **kwargs,
    )


class TestDownloadCompletesWhenTheArchiveDoes:
    def test_a_tool_that_never_exits_no_longer_hangs_the_run(
        self, tmp_path, hanging_datasets
    ):
        start = time.monotonic()
        _download(tmp_path, hanging_datasets, stall_timeout=60)
        elapsed = time.monotonic() - start

        assert elapsed < 30, f"waited {elapsed:.0f}s for a finished download"
        assert (tmp_path / "db" / "chuvi.fa").is_file()
        assert (tmp_path / "db" / "chuvi.csv").is_file()

    def test_the_downloaded_data_is_actually_used(self, tmp_path, hanging_datasets):
        _download(tmp_path, hanging_datasets, stall_timeout=60)
        assert "MKALLVGTSG" in (tmp_path / "db" / "chuvi.fa").read_text()

    def test_a_well_behaved_tool_is_unaffected(self, tmp_path, exiting_datasets):
        start = time.monotonic()
        _download(tmp_path, exiting_datasets, stall_timeout=60)
        assert time.monotonic() - start < 10
        assert (tmp_path / "db" / "chuvi.fa").is_file()


class TestArchiveCompleteness:
    def test_a_truncated_archive_is_not_accepted_as_success(self, tmp_path):
        """A half-written file must never be mistaken for a finished download."""
        script = tmp_path / "datasets"
        script.write_text(
            "#!/bin/sh\n"
            'while [ $# -gt 0 ]; do if [ "$1" = "--filename" ]; then out="$2"; fi;'
            " shift; done\n"
            'printf "PK\\003\\004 truncated" > "$out"\n'
            "sleep 300\n"
        )
        script.chmod(0o755)

        start = time.monotonic()
        with pytest.raises(RuntimeError, match="stalled"):
            _download(tmp_path, str(script), stall_timeout=5, attempts=1)
        # It gave up on the stall instead of accepting the partial file.
        assert time.monotonic() - start < 40


class TestOutputsAreCleanedUp:
    """Only the database and its provenance are kept, not the raw download."""

    def test_the_zip_and_extracted_directory_are_removed(
        self, tmp_path, exiting_datasets
    ):
        _download(tmp_path, exiting_datasets, stall_timeout=60)
        produced = sorted(p.name for p in (tmp_path / "db").iterdir())
        assert produced == ["chuvi.csv", "chuvi.fa", "chuvi.log", "chuvi.tracking.tsv"]

    @pytest.mark.skipif(
        not binaries_available("cd-hit"), reason="requires cd-hit on PATH"
    )
    def test_clustering_keeps_the_cluster_composition(self, tmp_path, exiting_datasets):
        _download(tmp_path, exiting_datasets, stall_timeout=60, cluster=True)
        clusters = tmp_path / "db" / "chuvi.clstr"
        assert clusters.is_file()
        assert ">Cluster" in clusters.read_text()

    def test_keep_download_preserves_them(self, tmp_path, exiting_datasets):
        _download(tmp_path, exiting_datasets, stall_timeout=60, keep_download=True)
        produced = sorted(p.name for p in (tmp_path / "db").iterdir())
        assert "chuvi.zip" in produced
        assert "chuvi_ncbi" in produced

    def test_the_tracking_table_covers_every_downloaded_accession(
        self, tmp_path, exiting_datasets
    ):
        _download(tmp_path, exiting_datasets, stall_timeout=60)
        table = pd.read_csv(tmp_path / "db" / "chuvi.tracking.tsv", sep="\t")
        assert list(table.Accession) == ["YP_000001.1"]
        assert set(table.Status) == {"kept"}


class TestHttp2Fallback:
    """NCBI resets HTTP/2 streams on large transfers; HTTP/1.1 completes them."""

    def test_a_stream_reset_is_retried_over_http1(self, tmp_path):
        package = tmp_path / "package.zip"
        _write_package(str(package))
        script = tmp_path / "datasets"
        script.write_text(
            "#!/bin/sh\n"
            'while [ $# -gt 0 ]; do if [ "$1" = "--filename" ]; then out="$2"; fi;'
            " shift; done\n"
            # Succeeds only once the client has been told to avoid HTTP/2.
            'if [ -z "$GODEBUG" ]; then\n'
            '  echo "Error: Download error: stream error: stream ID 3;'
            ' INTERNAL_ERROR; received from peer" >&2\n'
            "  exit 1\n"
            "fi\n"
            f'cp "{package}" "$out"\n'
        )
        script.chmod(0o755)

        _download(tmp_path, str(script), stall_timeout=30, attempts=3)
        assert (tmp_path / "db" / "chuvi.fa").is_file()


class TestBrokenArchives:
    """A download is only successful if what it left on disk can be read."""

    def test_a_clean_exit_with_a_broken_archive_is_not_accepted(self, tmp_path):
        """NCBI sometimes delivers a full-size package that is not a valid zip."""
        script = tmp_path / "datasets"
        script.write_text(
            "#!/bin/sh\n"
            'while [ $# -gt 0 ]; do if [ "$1" = "--filename" ]; then out="$2"; fi;'
            " shift; done\n"
            'printf "not a zip at all" > "$out"\n'
            'printf "Downloading: %s    118MB invalid zip archive\\n" "$out"\n'
            "exit 0\n"  # <- reports success anyway
        )
        script.chmod(0o755)

        with pytest.raises(RuntimeError) as failure:
            _download(tmp_path, str(script), stall_timeout=5, attempts=1)
        assert "archive" in str(failure.value)

    def test_a_leftover_archive_is_discarded_before_downloading(self, tmp_path):
        """An interrupted run must not leave something to download into."""
        outdir = tmp_path / "db"
        outdir.mkdir()
        stale = outdir / "chuvi.zip"
        stale.write_bytes(b"leftover from an interrupted run")

        package = tmp_path / "package.zip"
        _write_package(str(package))
        script = tmp_path / "datasets"
        # Refuses to run if the target already exists, which is what a resuming
        # client effectively does with a stale file.
        script.write_text(
            "#!/bin/sh\n"
            'while [ $# -gt 0 ]; do if [ "$1" = "--filename" ]; then out="$2"; fi;'
            " shift; done\n"
            'if [ -e "$out" ]; then echo "target exists" >&2; exit 1; fi\n'
            f'cp "{package}" "$out"\n'
        )
        script.chmod(0o755)

        _download(tmp_path, str(script), stall_timeout=30, attempts=1)
        assert (outdir / "chuvi.fa").is_file()
