"""Unit tests for eefinder.versions."""

from __future__ import annotations

import platform
import textwrap

import numpy
import pandas

from eefinder.versions import (
    DEPENDENCY_NAMES,
    _tool_version,
    collect_dependency_versions,
    collect_system_info,
    find_env_yml,
    parse_env_versions,
)

ENV_YML = textwrap.dedent("""\
    name: EEfinder
    dependencies:
      - bedtools=2.27.1
      - blast=2.5.0
      - diamond=2.0.15
      - python=3.9.0
      - pip:
        - numpy==1.23.1
        - pandas==1.4.2
    """)


def test_parse_env_versions_handles_conda_and_pip_pins(tmp_path):
    env = tmp_path / "env.yml"
    env.write_text(ENV_YML)

    assert parse_env_versions(env) == {
        "python": "3.9.0",
        "numpy": "1.23.1",
        "pandas": "1.4.2",
        "bedtools": "2.27.1",
        "blast": "2.5.0",
        "diamond": "2.0.15",
    }


def test_collect_flags_matches_and_mismatches(tmp_path):
    # Pin python/numpy to the actually-installed versions (-> ok) and pandas to
    # a bogus one (-> mismatch); the external tools are left unpinned here.
    env = tmp_path / "env.yml"
    env.write_text(textwrap.dedent(f"""\
            dependencies:
              - python={platform.python_version()}
              - pip:
                - numpy=={numpy.__version__}
                - pandas==0.0.0
            """))

    deps = {dep.name: dep for dep in collect_dependency_versions(env)}
    assert deps["python"].status == "ok"
    assert deps["numpy"].status == "ok"
    assert deps["pandas"].status == "mismatch"
    assert deps["pandas"].detected == pandas.__version__
    # Tools absent from this env.yml are either unpinned (present on PATH) or
    # not-found (absent) -- never a false "ok"/"mismatch".
    assert deps["bedtools"].status in {"unpinned", "not-found"}


def test_collect_prefix_pin_matches_patch_version(tmp_path):
    # A major.minor pin (e.g. "python=3.10") should match the running patch
    # release ("3.10.x"), not report a spurious mismatch.
    major_minor = ".".join(platform.python_version().split(".")[:2])
    env = tmp_path / "env.yml"
    env.write_text(f"dependencies:\n  - python={major_minor}\n")

    deps = {dep.name: dep for dep in collect_dependency_versions(env)}
    assert deps["python"].status == "ok"


def test_collect_reports_all_dependencies_in_order(tmp_path):
    names = [dep.name for dep in collect_dependency_versions(None)]
    assert names == list(DEPENDENCY_NAMES)


def test_collect_without_env_yml_is_unpinned(tmp_path):
    deps = {dep.name: dep for dep in collect_dependency_versions(None)}
    # python is always detectable but has no expected version to compare to.
    assert deps["python"].detected == platform.python_version()
    assert deps["python"].expected is None
    assert deps["python"].status == "unpinned"


def test_collect_system_info_is_populated():
    info = collect_system_info()
    # All fields are non-empty strings describing the host/OS.
    assert info.operating_system
    assert info.machine
    assert isinstance(info.hostname, str)
    assert isinstance(info.user, str)


def test_tool_version_missing_binary_returns_none():
    assert _tool_version("eefinder_no_such_binary_xyz --version") is None


def test_find_env_yml_honours_override(tmp_path, monkeypatch):
    env = tmp_path / "env.yml"
    env.write_text(ENV_YML)

    monkeypatch.setenv("EEFINDER_ENV_YML", str(env))
    assert find_env_yml() == env

    monkeypatch.setenv("EEFINDER_ENV_YML", str(tmp_path / "missing.yml"))
    assert find_env_yml() is None
