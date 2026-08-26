from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("eefinder")
except PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "unknown"
