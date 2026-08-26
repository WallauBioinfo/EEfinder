"""Resolve NCBI taxonomy so a subtree can be left out of a download.

The ``datasets`` CLI can download a taxon, but it cannot download a taxon
*minus* one of its branches: its ``virus genome taxon`` subcommand has no
exclusion flag. That matters because SARS-CoV-2 (tax id 2697049) is 61% of every
viral record in GenBank, so ``get-databases virus --all-sequences`` spends most
of its time and disk on one virus.

The way out is to never ask for the branch in the first place. For a root taxon
``R`` and an excluded taxon ``X`` below it, walk the lineage ``R -> ... -> X``
and, at each step, keep every child except the one that leads to ``X``. The
union of those siblings is exactly ``subtree(R) - subtree(X)``, and it is small
enough to hand back to ``datasets`` through ``--inputfile``: excluding
SARS-CoV-2 from all viruses yields 63 taxa.

One caveat is inherent to the method: records attached **directly** to an
internal node on that path (``"Betacoronavirus sp."`` and friends) belong to no
child, so they are not fetched. Measured against NCBI's own totals, that is 273
of 5.8 million records (0.005%) for the SARS-CoV-2 case -- all of them
unclassified entries whose genus and family EEfinder could not resolve anyway.
"""

from __future__ import annotations

import json
import shlex
import subprocess
from typing import NamedTuple

from eefinder.log import logger

#: The ``datasets`` executable.
DATASETS_BINARY = "datasets"

#: Taxa left out of a viral download unless the user asks otherwise.
#: 2697049 is SARS-CoV-2.
DEFAULT_VIRUS_EXCLUSIONS = ("2697049",)

#: Value of ``--exclude-taxon`` that turns the default exclusions off.
NO_EXCLUSION = "none"

#: Maximum number of taxa the ``datasets`` CLI accepts in an ``--inputfile``.
INPUTFILE_LIMIT = 100


class TaxonNode(NamedTuple):
    """One node of the NCBI taxonomy, as ``datasets`` reports it."""

    tax_id: int
    name: str
    #: Ancestors, root first (the lineage NCBI returns in ``parents``).
    parents: "tuple[int, ...]"
    #: Direct children only.
    children: "tuple[int, ...]"


class Expansion(NamedTuple):
    """The taxa to download, after excluded branches were pruned away."""

    #: Tax ids to hand to ``datasets``. A single-element tuple means no
    #: exclusion applied and the root can be requested directly.
    taxa: "tuple[int, ...]"
    root: TaxonNode
    #: Excluded taxa that were actually found below the root.
    excluded: "tuple[TaxonNode, ...]"
    #: Excluded taxa that are not below the root, so nothing was pruned.
    not_below_root: "tuple[TaxonNode, ...]"

    @property
    def pruned(self) -> bool:
        """Whether anything was actually excluded."""
        return bool(self.excluded)


def summarize_taxa(
    taxa: "list[str] | tuple[str, ...]",
    datasets_bin: str = DATASETS_BINARY,
) -> "dict[int, TaxonNode]":
    """Look up taxa by id or scientific name, returning them by tax id.

    Parameters
    ----------
    taxa : list[str]
        Tax ids or scientific names. Passed to ``datasets`` in one call.
    datasets_bin : str
        Path/name of the ``datasets`` executable.

    Returns
    -------
    dict[int, TaxonNode]
        One entry per taxon the CLI recognised.

    Raises
    ------
    RuntimeError
        If the ``datasets`` call fails.
    """
    if not taxa:
        return {}
    joined = ",".join(str(t) for t in taxa)
    command = (
        f"{datasets_bin} summary taxonomy taxon {shlex.quote(joined)} --as-json-lines"
    )
    logger.debug(f"taxonomy command: {command}")
    result = subprocess.run(
        shlex.split(command),
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"could not look up taxonomy for {joined!r}: "
            f"{(result.stderr or '').strip() or 'datasets exited ' + str(result.returncode)}"
        )
    nodes: "dict[int, TaxonNode]" = {}
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError:
            continue
        record = payload.get("taxonomy", payload)
        tax_id = record.get("tax_id")
        if tax_id is None:
            continue
        name = (record.get("current_scientific_name") or {}).get("name", "")
        nodes[int(tax_id)] = TaxonNode(
            tax_id=int(tax_id),
            name=name,
            parents=tuple(record.get("parents") or ()),
            children=tuple(record.get("children") or ()),
        )
    return nodes


class _Taxonomy:
    """Lazily fetched taxonomy nodes, cached across an expansion."""

    def __init__(self, datasets_bin: str = DATASETS_BINARY) -> None:
        self.datasets_bin = datasets_bin
        self._cache: "dict[int, TaxonNode]" = {}

    def node(self, tax_id: int) -> TaxonNode:
        """Return one node, fetching it unless it is already cached."""
        if tax_id not in self._cache:
            found = summarize_taxa([str(tax_id)], self.datasets_bin)
            if tax_id not in found:
                raise RuntimeError(f"NCBI taxonomy has no taxon {tax_id}")
            self._cache.update(found)
        return self._cache[tax_id]

    def prime(self, nodes: "dict[int, TaxonNode]") -> None:
        """Seed the cache with already-fetched nodes."""
        self._cache.update(nodes)


def resolve_taxon(
    taxon: str, datasets_bin: str = DATASETS_BINARY
) -> "TaxonNode | None":
    """Resolve a tax id or scientific name to a :class:`TaxonNode`."""
    found = summarize_taxa([str(taxon)], datasets_bin)
    if not found:
        return None
    key = None
    try:
        key = int(str(taxon))
    except ValueError:
        pass
    if key is not None and key in found:
        return found[key]
    return next(iter(found.values()))


def _siblings_along_path(
    taxonomy: "_Taxonomy", root_id: int, excluded_id: int
) -> "list[int]":
    """Taxa covering ``subtree(root) - subtree(excluded)``.

    Walks the lineage from ``root`` down to ``excluded`` and, at each node,
    collects every child except the one on the path.
    """
    excluded = taxonomy.node(excluded_id)
    lineage = list(excluded.parents)
    if root_id not in lineage:
        raise ValueError(f"{excluded_id} is not below {root_id}")
    path = lineage[lineage.index(root_id) :] + [excluded_id]

    keep: "list[int]" = []
    for parent_id, next_id in zip(path, path[1:]):
        parent = taxonomy.node(parent_id)
        kept = [child for child in parent.children if child != next_id]
        logger.debug(
            f"  {parent_id} ({parent.name}): {len(parent.children)} child(ren), "
            f"keeping {len(kept)}"
        )
        keep.extend(kept)
    return keep


def expand_taxon_excluding(
    root: str,
    excluded: "list[str] | tuple[str, ...]" = (),
    datasets_bin: str = DATASETS_BINARY,
) -> Expansion:
    """Expand ``root`` into the taxa to download, minus the excluded branches.

    Parameters
    ----------
    root : str
        Tax id or scientific name of the taxon being downloaded.
    excluded : list[str]
        Tax ids or scientific names to leave out. Entries that are not below
        ``root`` are reported in :attr:`Expansion.not_below_root` and otherwise
        ignored -- excluding something that was never going to be downloaded is
        not an error.
    datasets_bin : str
        Path/name of the ``datasets`` executable.

    Returns
    -------
    Expansion
        ``taxa`` is ``(root,)`` when nothing was pruned.

    Raises
    ------
    RuntimeError
        If ``root`` or an excluded taxon cannot be resolved, or if the root is
        itself excluded (which would leave nothing to download).
    """
    taxonomy = _Taxonomy(datasets_bin)
    root_node = resolve_taxon(root, datasets_bin)
    if root_node is None:
        raise RuntimeError(f"NCBI taxonomy has no taxon {root!r}")
    taxonomy.prime({root_node.tax_id: root_node})

    wanted = [str(t) for t in excluded if str(t).strip()]
    if not wanted:
        return Expansion((root_node.tax_id,), root_node, (), ())

    found = summarize_taxa(wanted, datasets_bin)
    if len(found) < len(set(wanted)):
        known = {node.name.lower() for node in found.values()} | {str(t) for t in found}
        missing = [t for t in wanted if t.lower() not in known and t not in known]
        if missing:
            raise RuntimeError(
                "NCBI taxonomy has no taxon " + ", ".join(repr(m) for m in missing)
            )
    taxonomy.prime(found)

    keep = [root_node.tax_id]
    applied: "list[TaxonNode]" = []
    outside: "list[TaxonNode]" = []
    for node in found.values():
        if node.tax_id == root_node.tax_id:
            raise RuntimeError(
                f"cannot exclude {node.name} ({node.tax_id}): it is the taxon "
                "being downloaded, so nothing would be left"
            )
        if root_node.tax_id not in node.parents:
            outside.append(node)
            continue
        logger.debug(f"Excluding {node.name} ({node.tax_id}) from the download")
        expanded: "list[int]" = []
        touched = False
        for current in keep:
            if current == node.tax_id:
                touched = True
                continue
            if current in node.parents:
                expanded.extend(_siblings_along_path(taxonomy, current, node.tax_id))
                touched = True
            else:
                expanded.append(current)
        keep = expanded
        if touched:
            applied.append(node)
        else:
            outside.append(node)

    # Deduplicate while keeping a stable order for reproducible commands.
    seen: "set[int]" = set()
    taxa = tuple(t for t in keep if not (t in seen or seen.add(t)))
    if not taxa:
        raise RuntimeError(
            "the exclusions leave no taxa to download; narrow -tx or drop one"
        )
    return Expansion(taxa, root_node, tuple(applied), tuple(outside))


def batch_taxa(
    taxa: "tuple[int, ...] | list[int]", limit: int = INPUTFILE_LIMIT
) -> "list[list[int]]":
    """Split ``taxa`` into groups the ``datasets`` ``--inputfile`` accepts."""
    ids = list(taxa)
    return [ids[i : i + limit] for i in range(0, len(ids), limit)] or [[]]
