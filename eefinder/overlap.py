"""Resolve elements tagged ``overlaped`` (the ``tag`` column) by a chosen strategy.

Elements are flagged ``overlaped`` by :class:`~eefinder.tag_elements.TagElements`
when they overlap (within a small margin, same contig) an element assigned to a
*different* family. This module decides what to do with them:

* ``keep``    -- keep every element (no filtering);
* ``longest`` -- among elements that physically overlap, keep the longest and
  drop the shorter ones;
* ``targets`` -- resolve each overlap *cluster* by a supplied family list, given
  as **either** a keep-list or a drop-list:

  - keep-list (``target_families``): a cluster that contains at least one
    target-family member keeps its target-family elements and drops the rest; a
    cluster with no target member is left untouched (every element kept);
  - drop-list (``non_target_families``): the listed families are dropped from
    each cluster, unless *every* member of the cluster is a listed family (in
    which case the cluster is left untouched, never wiped).

The filtered-out elements are not discarded: :class:`FilterOverlap` writes them
to a ``tmp_outputs/`` directory so they remain available for inspection.
"""

from __future__ import annotations

import re

import pandas as pd

#: Accepted values for the ``--overlap`` option.
OVERLAP_STRATEGIES = ("keep", "longest", "targets")

_OVERLAPED_TAG = "overlaped"


def _element_length(element_id: str) -> int:
    """Length (in bp) encoded in a ``contig:start-end`` element id."""
    _, _, coords = element_id.rpartition(":")
    start, _, end = coords.partition("-")
    return int(end) - int(start)


def _overlap_clusters(overlaped: pd.DataFrame) -> list[set[str]]:
    """Group overlaping elements into connected clusters.

    Two elements belong to the same cluster when one lists the other in its
    ``Overlaped_Element_ID`` (transitively): the clusters are the connected
    components of the overlap graph, restricted to elements tagged ``overlaped``.

    Parameters
    ----------
    overlaped : pandas.DataFrame
        The ``overlaped`` rows of a taxonomy table (``Element-ID`` and
        ``Overlaped_Element_ID`` columns).

    Returns
    -------
    list[set[str]]
        One set of ``Element-ID``s per connected cluster.
    """
    ids = set(overlaped["Element-ID"])
    adjacency: dict[str, set[str]] = {eid: set() for eid in ids}
    for _, row in overlaped.iterrows():
        eid = row["Element-ID"]
        for partner in str(row["Overlaped_Element_ID"]).split(","):
            if partner in ids:
                adjacency[eid].add(partner)
                adjacency[partner].add(eid)

    clusters: list[set[str]] = []
    seen: set[str] = set()
    for start in ids:
        if start in seen:
            continue
        stack = [start]
        cluster: set[str] = set()
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            cluster.add(node)
            stack.extend(adjacency[node] - seen)
        clusters.append(cluster)
    return clusters


def elements_to_remove(
    df: pd.DataFrame,
    strategy: str,
    target_families: list[str],
    non_target_families: "list[str] | None" = None,
) -> set[str]:
    """Return the ``Element-ID``s to drop for the given overlap strategy.

    Only elements tagged ``overlaped`` are ever considered for removal; unique
    elements are always kept.

    Parameters
    ----------
    df : pandas.DataFrame
        Taxonomy table (must have ``Element-ID``, ``Family``, ``tag`` and
        ``Overlaped_Element_ID`` columns).
    strategy : str
        One of :data:`OVERLAP_STRATEGIES`.
    target_families : list[str]
        Families to *keep*, used only when ``strategy == "targets"``. Mutually
        exclusive with ``non_target_families``.
    non_target_families : list[str], optional
        Families to *drop*, used only when ``strategy == "targets"`` and
        ``target_families`` is empty. Mutually exclusive with
        ``target_families``.

    Returns
    -------
    set[str]
        The ``Element-ID``s to remove.
    """
    if strategy == "keep":
        return set()

    overlaped = df[df["tag"] == _OVERLAPED_TAG]

    if strategy == "targets":
        targets = set(target_families)
        non_targets = set(non_target_families or [])
        family_by_id = dict(zip(df["Element-ID"], df["Family"]))
        removed: set[str] = set()
        for cluster in _overlap_clusters(overlaped):
            families = {family_by_id.get(eid) for eid in cluster}
            if targets:
                # Keep-list mode: resolve only clusters that contain a target
                # family, dropping every non-target member; a cluster with no
                # target member is left untouched (all kept).
                if not (families & targets):
                    continue
                removed.update(
                    eid for eid in cluster if family_by_id.get(eid) not in targets
                )
            else:
                # Drop-list mode: drop the listed families, but never wipe a
                # cluster whose every member is a non-target family.
                if families <= non_targets:
                    continue
                removed.update(
                    eid for eid in cluster if family_by_id.get(eid) in non_targets
                )
        return removed

    if strategy == "longest":
        length_by_id = {eid: _element_length(eid) for eid in df["Element-ID"]}
        removed: set[str] = set()
        for _, row in overlaped.iterrows():
            own_length = length_by_id[row["Element-ID"]]
            partners = [p for p in str(row["Overlaped_Element_ID"]).split(",") if p]
            # Drop this element if any element it overlaps is strictly longer.
            if any(length_by_id.get(p, 0) > own_length for p in partners):
                removed.add(row["Element-ID"])
        return removed

    raise ValueError(f"Unknown overlap strategy: {strategy!r}")


class FilterOverlap:
    """Filter overlaping elements out of an ``(fasta, taxonomy)`` result pair.

    The kept elements are written back over ``fasta_file`` and ``tax_file``; the
    removed elements are written to ``removed_fasta`` / ``removed_tax``. Runs on
    instantiation.

    Parameters
    ----------
    fasta_file : str
        EE FASTA to filter in place (headers are ``{prefix}/{Element-ID}``).
    tax_file : str
        Taxonomy TSV to filter in place.
    strategy : str
        One of :data:`OVERLAP_STRATEGIES`.
    target_families : list[str]
        Families to keep for the ``targets`` strategy (keep-list).
    removed_fasta : str
        Path to write the filtered-out sequences.
    removed_tax : str
        Path to write the filtered-out taxonomy rows.
    non_target_families : list[str], optional
        Families to drop for the ``targets`` strategy (drop-list), used when
        ``target_families`` is empty.
    """

    def __init__(
        self,
        fasta_file: str,
        tax_file: str,
        strategy: str,
        target_families: list[str],
        removed_fasta: str,
        removed_tax: str,
        non_target_families: "list[str] | None" = None,
    ) -> None:
        self.fasta_file = fasta_file
        self.tax_file = tax_file
        self.strategy = strategy
        self.target_families = target_families
        self.non_target_families = non_target_families or []
        self.removed_fasta = removed_fasta
        self.removed_tax = removed_tax

        self.filter_overlap()

    def filter_overlap(self) -> None:
        """Split the taxonomy table and FASTA into kept vs removed elements."""
        df = pd.read_csv(self.tax_file, sep="\t")
        removed_ids = elements_to_remove(
            df, self.strategy, self.target_families, self.non_target_families
        )

        kept = df[~df["Element-ID"].isin(removed_ids)]
        removed = df[df["Element-ID"].isin(removed_ids)]
        kept.to_csv(self.tax_file, sep="\t", index=False)
        removed.to_csv(self.removed_tax, sep="\t", index=False)

        self._split_fasta(removed_ids)

    def _split_fasta(self, removed_ids: set[str]) -> None:
        """Route each FASTA record to the kept or removed file (verbatim)."""
        kept_lines: list[str] = []
        removed_lines: list[str] = []
        bucket = kept_lines
        with open(self.fasta_file) as fasta_in:
            for line in fasta_in:
                if line.startswith(">"):
                    # Strip the "{prefix}/" so the header matches the Element-ID.
                    element_id = re.sub(r".*/", "", line[1:].rstrip("\n"))
                    bucket = removed_lines if element_id in removed_ids else kept_lines
                bucket.append(line)

        with open(self.fasta_file, "w") as kept_out:
            kept_out.writelines(kept_lines)
        with open(self.removed_fasta, "w") as removed_out:
            removed_out.writelines(removed_lines)
