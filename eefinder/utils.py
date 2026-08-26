"""Filesystem and run-bookkeeping helpers for the EEfinder pipeline.

The timing/summary structures emitted at the end of a run are modelled as
dataclasses (:class:`StepInfo`, :class:`RunArguments`, :class:`RunInfo`) so the
run log has a single, typed source of truth. Use :func:`dataclasses.asdict` to
serialise them to JSON.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from eefinder.versions import DependencyVersion, SystemInfo

#: Timestamp format used throughout the run log.
TIME_FORMAT = "%Y-%m-%d %H:%M:%S"


def check_outdir(outdir: str) -> str:
    """Create the output directory and return its normalised path.

    Parameters
    ----------
    outdir : str
        Desired output directory, with or without a trailing slash.

    Returns
    -------
    str
        The directory path without a trailing slash. The directory (and any
        missing parents) is created if it does not already exist.

    Example
    -------
    >>> check_outdir("results/")  # doctest: +SKIP
    'results'
    """
    outdir = re.sub("/$", "", outdir)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    return outdir


def _format_timestamp(epoch: float) -> str:
    """Render a POSIX timestamp using :data:`TIME_FORMAT`."""
    return datetime.fromtimestamp(epoch).strftime(TIME_FORMAT)


def _elapsed_minutes(start_time: float, end_time: float) -> str:
    """Return the elapsed wall-clock time in minutes, formatted to 4 decimals."""
    return f"{(end_time - start_time) / 60:.4f}"


@dataclass
class StepInfo:
    """Timing and description record for a single pipeline step."""

    step: str
    start_time: str
    end_time: str
    total_time_minutes: str
    message: str

    @classmethod
    def from_times(
        cls, step: str, start_time: float, end_time: float, message: str
    ) -> "StepInfo":
        """Build a :class:`StepInfo` from POSIX start/end timestamps.

        Parameters
        ----------
        step : str
            Human-readable name of the pipeline step.
        start_time, end_time : float
            POSIX timestamps (e.g. from :func:`time.time`) bounding the step.
        message : str
            Description of what the step produced.

        Returns
        -------
        StepInfo
        """
        return cls(
            step=step,
            start_time=_format_timestamp(start_time),
            end_time=_format_timestamp(end_time),
            total_time_minutes=_elapsed_minutes(start_time, end_time),
            message=message,
        )


@dataclass
class RunArguments:
    """The full set of resolved command-line arguments for a run."""

    genome_file: str
    prefix: str
    outdir: str
    database: str
    dbmetadata: str
    baits: str
    mode: str
    length: int
    flank: int
    limit: int
    range_junction: int
    mask_per: int
    clean_masked: bool
    threads: int
    removetmp: bool
    index_databases: bool
    merge_level: str
    analysis: str
    overlap: str
    target_families: list
    non_target_families: list
    translation_method: str


@dataclass
class DownloadArguments:
    """The resolved arguments of a ``get-databases`` run."""

    dataset: str
    taxon: str
    outdir: str
    prefix: str
    refseq: bool
    exclude_uninformative: bool
    standardize_proteins: bool
    cluster: bool
    #: Download robustness settings (retries and stall detection).
    attempts: int = 1
    stall_timeout: float = 0.0
    keep_download: bool = False
    released_before: str = ""
    #: Taxa left out of the download entirely (never requested from NCBI).
    exclude_taxa: str = ""


@dataclass
class SequenceCounts:
    """How many sequences a ``get-databases`` run kept vs dropped."""

    downloaded: int
    excluded_uninformative: int
    clustered_identical: int
    dropped_standardization: int
    kept: int


@dataclass
class DownloadInfo:
    """Top-level ``get-databases`` summary serialised to ``{prefix}.log``."""

    eefinder_version: str
    arguments: DownloadArguments
    sequence_counts: SequenceCounts
    start_time: str
    end_time: str
    total_time_minutes: str
    steps_information: list[StepInfo] = field(default_factory=list)

    @classmethod
    def from_run(
        cls,
        eefinder_version: str,
        arguments: DownloadArguments,
        sequence_counts: SequenceCounts,
        start_time: float,
        end_time: float,
        steps_information: list[StepInfo],
    ) -> "DownloadInfo":
        """Assemble the download summary from metadata, timestamps and steps.

        Parameters
        ----------
        eefinder_version : str
            Version of EEfinder that produced the download.
        arguments : DownloadArguments
            Resolved arguments used for the download.
        sequence_counts : SequenceCounts
            Kept vs dropped sequence tallies.
        start_time, end_time : float
            POSIX timestamps bounding the whole download.
        steps_information : list[StepInfo]
            One :class:`StepInfo` per download phase.

        Returns
        -------
        DownloadInfo
        """
        return cls(
            eefinder_version=eefinder_version,
            arguments=arguments,
            sequence_counts=sequence_counts,
            start_time=_format_timestamp(start_time),
            end_time=_format_timestamp(end_time),
            total_time_minutes=_elapsed_minutes(start_time, end_time),
            steps_information=steps_information,
        )


@dataclass
class RunInfo:
    """Top-level run summary serialised to ``eefinder.log``."""

    eefinder_version: str
    system: "SystemInfo"
    arguments: RunArguments
    dependencies: list["DependencyVersion"]
    start_time: str
    end_time: str
    total_time_minutes: str
    steps_information: list[StepInfo] = field(default_factory=list)

    @classmethod
    def from_run(
        cls,
        eefinder_version: str,
        system: "SystemInfo",
        arguments: RunArguments,
        dependencies: list["DependencyVersion"],
        start_time: float,
        end_time: float,
        steps_information: list[StepInfo],
    ) -> "RunInfo":
        """Assemble the run summary from run metadata, timestamps and step info.

        Parameters
        ----------
        eefinder_version : str
            Version of EEfinder that produced the run.
        system : SystemInfo
            Operating system and host context of the run.
        arguments : RunArguments
            Resolved arguments used for the run.
        dependencies : list[DependencyVersion]
            Detected vs env.yml-pinned versions of the dependencies.
        start_time, end_time : float
            POSIX timestamps bounding the whole run.
        steps_information : list[StepInfo]
            One :class:`StepInfo` per executed pipeline step.

        Returns
        -------
        RunInfo
        """
        return cls(
            eefinder_version=eefinder_version,
            system=system,
            arguments=arguments,
            dependencies=dependencies,
            start_time=_format_timestamp(start_time),
            end_time=_format_timestamp(end_time),
            total_time_minutes=_elapsed_minutes(start_time, end_time),
            steps_information=steps_information,
        )
