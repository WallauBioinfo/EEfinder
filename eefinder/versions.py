"""Detect runtime dependency versions and compare them against env.yml pins.

Used to enrich ``eefinder.log`` with the versions of the external tools
(bedtools, BLAST, DIAMOND) and Python libraries (python, numpy, pandas) that a
run actually used, flagging any that differ from the versions pinned in
``env.yml`` (i.e. the versions EEfinder was validated with).
"""

from __future__ import annotations

import getpass
import os
import platform
import re
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy
import pandas

#: External binaries and the command used to query each one's version.
_TOOL_COMMANDS = {
    "bedtools": "bedtools --version",
    "blast": "blastx -version",
    "diamond": "diamond version",
    "ncbi-datasets-cli": "datasets --version",
}

#: Order in which dependencies are reported in the run log.
DEPENDENCY_NAMES = (
    "python",
    "numpy",
    "pandas",
    "bedtools",
    "blast",
    "diamond",
    "ncbi-datasets-cli",
)

_VERSION_RE = re.compile(r"(\d+\.\d+(?:\.\d+)?)")

# Status values for a dependency-version comparison.
STATUS_OK = "ok"
STATUS_MISMATCH = "mismatch"
STATUS_NOT_FOUND = "not-found"
STATUS_UNPINNED = "unpinned"


@dataclass
class DependencyVersion:
    """Detected vs env.yml-declared version of a single dependency.

    Attributes
    ----------
    name : str
        Dependency name.
    detected : str | None
        Version found at runtime, or ``None`` if the tool is not on ``PATH``.
    expected : str | None
        Version pinned in ``env.yml``, or ``None`` if it is not pinned/known.
    status : str
        One of ``"ok"``, ``"mismatch"``, ``"not-found"`` or ``"unpinned"``.
    """

    name: str
    detected: str | None
    expected: str | None
    status: str


@dataclass
class SystemInfo:
    """Host/run context (operating system, machine, host and user)."""

    operating_system: str
    machine: str
    hostname: str
    user: str


def collect_system_info() -> SystemInfo:
    """Gather the operating system and host context for the run log."""
    try:
        user = getpass.getuser()
    except Exception:  # pragma: no cover - environment-dependent
        user = "unknown"
    return SystemInfo(
        operating_system=platform.platform(),
        machine=platform.machine(),
        hostname=platform.node(),
        user=user,
    )


def find_env_yml() -> Path | None:
    """Locate ``env.yml``.

    Uses the ``EEFINDER_ENV_YML`` environment variable if set, otherwise looks
    for the file shipped alongside the package (available when running from a
    cloned/editable checkout). Returns ``None`` when it cannot be found.
    """
    override = os.environ.get("EEFINDER_ENV_YML")
    if override:
        path = Path(override)
        return path if path.is_file() else None
    candidate = Path(__file__).resolve().parent.parent / "env.yml"
    return candidate if candidate.is_file() else None


def _tool_version(command: str) -> str | None:
    """Return the version string reported by an external tool, or ``None``."""
    try:
        result = subprocess.run(shlex.split(command), capture_output=True, text=True)
    except (FileNotFoundError, OSError):
        return None
    match = _VERSION_RE.search(f"{result.stdout}\n{result.stderr}")
    return match.group(1) if match else None


def _detected_versions() -> dict[str, str | None]:
    """Detect the runtime version of every tracked dependency."""
    versions: dict[str, str | None] = {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
    }
    for name, command in _TOOL_COMMANDS.items():
        versions[name] = _tool_version(command)
    return versions


def parse_env_versions(env_yml: Path, names=DEPENDENCY_NAMES) -> dict[str, str]:
    """Parse the pinned versions of ``names`` from an ``env.yml`` file.

    Handles both conda (``- name=version``) and pip (``- name==version``) pins.
    """
    text = env_yml.read_text()
    pins: dict[str, str] = {}
    for name in names:
        match = re.search(rf"(?m)^\s*-\s*{re.escape(name)}\s*={{1,2}}\s*([\w.]+)", text)
        if match:
            pins[name] = match.group(1)
    return pins


def _normalize(version: str) -> str:
    """Strip cosmetic decorations so versions compare cleanly."""
    return version.strip().lstrip("v").rstrip("+")


def _versions_match(detected: str, expected: str) -> bool:
    """Whether a detected version satisfies an env.yml pin.

    A shorter pin is treated as a version prefix, so ``3.9`` matches ``3.9.23``
    while ``2.16.0`` does not match ``2.17.0``.
    """
    detected, expected = _normalize(detected), _normalize(expected)
    return detected == expected or detected.startswith(f"{expected}.")


def _status(detected: str | None, expected: str | None) -> str:
    if detected is None:
        return STATUS_NOT_FOUND
    if expected is None:
        return STATUS_UNPINNED
    return STATUS_OK if _versions_match(detected, expected) else STATUS_MISMATCH


def collect_dependency_versions(
    env_yml: Path | None = None,
) -> list[DependencyVersion]:
    """Return the detected vs expected version of every tracked dependency.

    Parameters
    ----------
    env_yml : Path | None
        Path to the ``env.yml`` providing the expected versions. When ``None``
        or missing, every dependency is reported with status ``"unpinned"``.
    """
    detected = _detected_versions()
    expected = parse_env_versions(env_yml) if env_yml and env_yml.is_file() else {}
    return [
        DependencyVersion(
            name=name,
            detected=detected.get(name),
            expected=expected.get(name),
            status=_status(detected.get(name), expected.get(name)),
        )
        for name in DEPENDENCY_NAMES
    ]
