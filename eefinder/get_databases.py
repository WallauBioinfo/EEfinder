"""Download EEfinder similarity-search databases via the NCBI datasets CLI.

EEfinder needs, as inputs to the ``screening`` command:

* a **protein database** FASTA (``-db``) plus a **metadata CSV** (``-mt``) with the
  columns ``Accession,Species,Genus,Family,Molecule_type,Protein,Host``; and
* a **host-gene baits** FASTA (``-bt``).

This module automates acquiring them from NCBI RefSeq — the manual procedure
described in the wiki (https://github.com/WallauBioinfo/EEfinder/wiki) — using
the `NCBI datasets <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/>`_ command
line tool (``ncbi-datasets-cli`` in ``env.yml``). Three dataset types are
supported:

* ``virus``    -- RefSeq viral proteins + metadata (the ``-db``/``-mt`` inputs);
* ``bacteria`` -- RefSeq bacterial proteins + metadata (bacteria screening mode);
* ``host``     -- RefSeq host proteins used as ``-bt`` baits (no metadata CSV).

The metadata CSV is rebuilt from two sources bundled in the datasets download:
the ``protein.faa`` headers (``Accession``, ``Protein`` product and ``Species``)
and the ``data_report.jsonl`` taxonomy report (``Genus`` and ``Family`` inferred
from the unranked lineage by ICTV name suffix, plus ``Host``), joined by
organism name. ``Molecule_type`` is not in the datasets report, so it is filled
from the bundled ICTV genome-composition table (``data/``) keyed by family.
"""

from __future__ import annotations

import glob
import json
import os
import re
import shlex
import shutil
import subprocess
import time
import zipfile
from datetime import datetime
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import NamedTuple

import pandas as pd
from Bio import SeqIO

from eefinder import __version__
from eefinder.log import logger
from eefinder.utils import (
    check_outdir,
    DownloadArguments,
    DownloadInfo,
    SequenceCounts,
    StepInfo,
)
from eefinder.normalization import standardize_protein, strip_bracket_tags
from eefinder.progress import (
    DEFAULT_ATTEMPTS,
    DEFAULT_STALL_TIMEOUT,
    STALLED,
    progress_bar,
    run_with_retries,
)
from eefinder.translation import cluster_proteins, parse_cdhit_clusters
from eefinder.taxon_exclusion import batch_taxa, expand_taxon_excluding

#: cd-hit executable used to collapse identical proteins (from ``cd-hit``).
CDHIT_BINARY = "cd-hit"

#: Dataset types this module can download.
DATASET_CHOICES = ("virus", "bacteria", "host")

#: Dataset types that also produce a metadata CSV.
_METADATA_DATASETS = ("virus", "bacteria")

#: Default NCBI taxon per dataset when ``--taxon`` is omitted. ``virus`` and
#: ``bacteria`` default to their whole-database roots (10239 = Viruses,
#: 2 = Bacteria); ``host`` has no default.
DEFAULT_TAXA = {"virus": "10239", "bacteria": "2"}

#: Product substrings whose proteins carry no taxonomic signal; optionally
#: dropped from a download via ``--exclude-uninformative``.
UNINFORMATIVE_PRODUCTS = ("hypothetical protein", "uncharacterized protein")

#: Columns of ``{prefix}.tracking.tsv``: the fate of every downloaded accession.
#: ``Organism_release_date`` is the earliest release date among the organism's
#: genome records -- the finest granularity available, since a protein in the
#: download cannot be traced back to the record it came from.
TRACKING_COLUMNS = [
    "Accession",
    "Species",
    "Protein_downloaded",
    "Protein_final",
    "Name_changed",
    "Status",
    "Reason",
    "Organism_release_date",
    "Cluster",
    "Cluster_representative",
]

#: Values of the tracking ``Status`` column.
STATUS_KEPT = "kept"
STATUS_REMOVED = "removed"

#: Values of the tracking ``Reason`` column, one per way a sequence can leave.
REASON_UNINFORMATIVE = "uninformative_product"
REASON_DUPLICATE = "identical_duplicate"
REASON_UNKNOWN_PRODUCT = "product_standardized_to_unknown"
REASON_ABSENT = "absent_from_final_database"
REASON_TOO_RECENT = "released_after_cutoff"

#: ``cd-hit`` silently throws away sequences of at most ``-l`` residues
#: (default 10), which is the usual reason a record vanishes without any step
#: reporting it. :func:`reconcile_tracking` catches that and anything like it.
CDHIT_MIN_LENGTH = 10

#: Columns of the metadata CSV consumed by the ``screening`` command (``-mt``).
METADATA_COLUMNS = [
    "Accession",
    "Species",
    "Genus",
    "Family",
    "Molecule_type",
    "Protein",
    "Host",
]

#: NCBI datasets CLI binary (from ``ncbi-datasets-cli``).
DATASETS_BINARY = "datasets"

#: Messages that mean NCBI reset the HTTP/2 stream mid-transfer. The download
#: itself is fine over HTTP/1.1, which the Go client can be told to use.
HTTP2_ERROR_MARKERS = ("stream error", "internal_error", "http2")

#: Environment that makes a Go binary fall back to HTTP/1.1.
HTTP1_ENV = {"GODEBUG": "http2client=0"}

#: Bundled ICTV family -> genome-composition table (used for ``Molecule_type``,
#: which the NCBI datasets report does not provide). Sourced from
#: https://ictv.global/virus-properties.
_ICTV_GENOME_TABLE = (
    Path(__file__).resolve().parent / "data" / ("ictv_genome_composition.tsv")
)


def _load_genome_composition() -> "dict[str, str]":
    """Load the ICTV ``family -> genome composition`` map from the data file."""
    table: dict[str, str] = {}
    if _ICTV_GENOME_TABLE.is_file():
        with open(_ICTV_GENOME_TABLE) as handle:
            next(handle, None)  # skip the header row
            for line in handle:
                family, _, genome = line.rstrip("\n").partition("\t")
                if family:
                    table[family] = genome
    return table


#: ICTV genome composition keyed by virus family.
GENOME_COMPOSITION = _load_genome_composition()


def molecule_type_for_family(family: str) -> str:
    """Return the ICTV genome composition (``Molecule_type``) for a family.

    Empty string when the family is unknown (e.g. an unclassified virus or any
    bacterial family, which the ICTV virus table does not cover).
    """
    return GENOME_COMPOSITION.get(family, "")


@dataclass
class ProteinHeader:
    """The three fields parsed out of a protein FASTA header."""

    accession: str
    product: str
    organism: str


@dataclass
class TaxonomyRecord:
    """Per-organism taxonomy pulled from a datasets ``data_report.jsonl``."""

    species: str
    genus: str
    family: str
    mol_type: str
    host: str


def parse_protein_header(header: str) -> ProteinHeader:
    """Parse a protein FASTA header into accession/product/organism.

    Handles both header styles EEfinder databases come in:

    * NCBI *datasets* CDS proteins ---
      ``YP_013613119.1:1-301 P1 [organism=Paris potyvirus 4] [isolate=YLJ]``:
      the organism is the ``[organism=...]`` tag and the product is the text
      before the first bracket group.
    * NCBI Virus RefSeq proteins ---
      ``YP_009664712.1 N protein [Bas-Congo tibrovirus]``: the organism is the
      trailing bare ``[...]`` group.

    In both cases the accession is the first whitespace-delimited token (kept
    whole, including any ``:start-end`` suffix, so it matches the FASTA id
    reported by BLAST).

    Parameters
    ----------
    header : str
        A FASTA header line, with or without a leading ``>``.

    Returns
    -------
    ProteinHeader
        ``organism``/``product`` are empty strings when the header lacks them.
    """
    header = header.lstrip(">").strip()
    accession, _, remainder = header.partition(" ")
    remainder = remainder.strip()

    organism = ""
    # The organism value may contain one level of nested "[...]" (a strain tag),
    # e.g. "[organism=Maize streak virus - A[South Africa]]".
    match = re.search(r"\[organism=((?:[^\[\]]|\[[^\[\]]*\])*)\]", remainder)
    if match:
        # datasets CDS format: product is everything before the first bracket.
        organism = match.group(1).strip()
        first_bracket = remainder.find("[")
        product = remainder[:first_bracket].strip()
    else:
        # RefSeq format: organism is the trailing bare "[...]" group, if any.
        trailing = re.search(r"\[([^\[\]]*)\]\s*$", remainder)
        if trailing:
            organism = trailing.group(1).strip()
            product = remainder[: trailing.start()].strip()
        else:
            product = remainder
    # Defensively drop any leaked "[key=value]" tag (e.g. "[organism=...]") the
    # branch above may have left behind for unusual header layouts.
    product = re.sub(r"\s+", " ", strip_bracket_tags(product)).strip()
    return ProteinHeader(accession=accession, product=product, organism=organism)


def _genus_family_from_lineage(lineage: "list[dict]") -> "tuple[str, str]":
    """Infer ``(genus, family)`` from an unranked datasets ``lineage``.

    The datasets ``data_report.jsonl`` lists the taxonomic lineage as
    ``{"name", "taxId"}`` entries **without** rank information, so ranks are
    inferred from the ICTV virus naming suffixes: families end in ``-viridae``
    and a genus is a single-word name ending in ``-virus``. The most specific
    (last) match of each is used.

    Returns empty strings for names that do not follow the ICTV conventions
    (e.g. unclassified viruses or phages).
    """
    genus = family = ""
    for node in lineage:
        name = node.get("name", "")
        lower = name.lower()
        if lower.endswith("viridae"):
            family = name
        elif " " not in name and lower.endswith("virus"):
            genus = name
    return genus, family


def parse_release_dates(report_path: str) -> "dict[str, str]":
    """Map each organism to the **earliest** release date of its records.

    A virus is represented by one record per genome (per segment, for segmented
    families), and those can be released years apart -- La Crosse virus has
    segments from 2002 and 2023. The earliest date is when the organism first
    appeared in RefSeq, which is what a cutoff is asking about.

    Parameters
    ----------
    report_path : str
        ``data_report.jsonl`` from the download.

    Returns
    -------
    dict[str, str]
        ``organism -> YYYY-MM-DD``; organisms without a date are absent.
    """
    dates: dict[str, str] = {}
    with open(report_path) as report:
        for line in report:
            line = line.strip()
            if not line:
                continue
            data = json.loads(line)
            organism = (data.get("virus") or data.get("organism") or {}).get(
                "organismName", ""
            )
            released = (data.get("releaseDate") or "")[:10]
            if not organism or not released:
                continue
            if organism not in dates or released < dates[organism]:
                dates[organism] = released
    return dates


def organisms_released_after(report_path: str, cutoff: str) -> "set[str]":
    """Organisms whose earliest record was released after ``cutoff``.

    These are the ones a dated build must leave out. The granularity is the
    organism, not the individual genome record: the protein FASTA carries no
    link back to the record a protein came from (its ``proteinCount`` ordering
    does not hold), so an organism that already existed at the cutoff is kept
    with all of its proteins.
    """
    return {
        organism
        for organism, released in parse_release_dates(report_path).items()
        if released > cutoff
    }


def parse_taxonomy_report(report_path: str) -> dict[str, TaxonomyRecord]:
    """Build an ``organism -> TaxonomyRecord`` map from a ``data_report.jsonl``.

    The datasets virus data report is a JSON-lines file with one record per
    genome. ``Species`` is the ``virus.organismName``, ``Genus``/``Family`` are
    inferred from the (unranked) ``virus.lineage`` via
    :func:`_genus_family_from_lineage`, and ``Host`` is the top-level
    ``host.organismName``. ``Molecule_type`` is left empty: the datasets report
    does not carry it. Fields are read defensively; the first record wins per
    organism.

    Parameters
    ----------
    report_path : str
        Path to the ``data_report.jsonl`` extracted from the datasets download.

    Returns
    -------
    dict[str, TaxonomyRecord]
        Keyed by organism (species) name.
    """
    records: dict[str, TaxonomyRecord] = {}
    with open(report_path) as report:
        for line in report:
            line = line.strip()
            if not line:
                continue
            data = json.loads(line)
            virus = data.get("virus") or data.get("organism") or {}
            organism = virus.get("organismName", "")
            if not organism or organism in records:
                continue

            genus, family = _genus_family_from_lineage(virus.get("lineage", []))
            host = (data.get("host") or {}).get("organismName", "")
            records[organism] = TaxonomyRecord(
                species=organism,
                genus=genus,
                family=family,
                mol_type="",  # not present in the datasets report
                host=host,
            )
    return records


def build_metadata_frame(
    protein_fasta: str,
    taxonomy: dict[str, TaxonomyRecord],
    standardize: bool = False,
    dataset: str = "virus",
    tracking: "dict[str, dict] | None" = None,
) -> pd.DataFrame:
    """Join protein headers with taxonomy into the screening metadata table.

    Parameters
    ----------
    protein_fasta : str
        Concatenated protein FASTA (its headers supply Accession/Protein/Species).
    taxonomy : dict[str, TaxonomyRecord]
        Organism -> taxonomy map from :func:`parse_taxonomy_report`.
    standardize : bool
        When ``True``, rewrite the ``Protein`` column with
        :func:`standardize_protein` (canonical names via the target-specific
        logic, special characters removed, leading letter capitalised). Records
        that standardise to ``"Unknown"`` (bare ``CDS``/``ORF`` directives) are
        **dropped** from the table.
    dataset : str
        The database target (``"virus"``/``"bacteria"``), selecting which
        standardisation logic :func:`standardize_protein` applies.
    tracking : dict[str, dict], optional
        Updated with the standardised name of each record and whether that
        renamed it (``RNA-dependent RNA polymerase`` -> ``RdRp``), plus the
        records dropped because they standardise to ``"Unknown"``.

    Returns
    -------
    pandas.DataFrame
        One row per (kept) protein, with :data:`METADATA_COLUMNS`.
    """
    rows = []
    with progress_bar(
        SeqIO.parse(protein_fasta, "fasta"),
        "Building metadata",
        length=count_fasta_records(protein_fasta),
    ) as records:
        for record in records:
            header = parse_protein_header(record.description)
            tax = taxonomy.get(
                header.organism,
                TaxonomyRecord(header.organism, "", "", "", ""),
            )
            # Molecule_type is absent from the datasets report, so it is looked
            # up from the ICTV genome-composition table by family.
            mol_type = molecule_type_for_family(tax.family)
            if standardize:
                protein = standardize_protein(header.product, mol_type, target=dataset)
                if tracking is not None and header.accession in tracking:
                    row = tracking[header.accession]
                    row["Protein_final"] = protein
                    row["Name_changed"] = "yes" if protein != header.product else "no"
                    if protein == "Unknown":
                        row["Status"] = STATUS_REMOVED
                        row["Reason"] = REASON_UNKNOWN_PRODUCT
                if protein == "Unknown":
                    continue  # bare CDS/ORF: drop the record entirely
            else:
                protein = header.product
                if tracking is not None and header.accession in tracking:
                    tracking[header.accession]["Protein_final"] = protein
                    tracking[header.accession]["Name_changed"] = "no"
            rows.append(
                {
                    "Accession": header.accession,
                    "Species": header.organism,
                    "Genus": tax.genus,
                    "Family": tax.family,
                    "Molecule_type": mol_type,
                    "Protein": protein,
                    "Host": tax.host,
                }
            )
    return pd.DataFrame(rows, columns=METADATA_COLUMNS)


def build_download_command(
    dataset: str,
    taxon: str,
    zip_path: str,
    refseq: bool = True,
    datasets_bin: str = DATASETS_BINARY,
    released_before: "str | None" = None,
    inputfile: "str | None" = None,
) -> str:
    """Compose the ``datasets download`` command for a dataset type.

    Parameters
    ----------
    dataset : str
        One of :data:`DATASET_CHOICES`.
    taxon : str
        NCBI taxon name or tax id to download (e.g. ``Flaviviridae``).
    zip_path : str
        Destination path for the downloaded ``ncbi_dataset`` zip.
    refseq : bool
        Restrict to RefSeq sequences (recommended).
    datasets_bin : str
        Path/name of the ``datasets`` executable.
    released_before : str, optional
        Only include records released on or before this date (``YYYY-MM-DD``).
        The ``datasets`` CLI implements this for genome downloads
        (bacteria/host); its **virus** subcommand has no such flag, so for
        viruses the cutoff is applied afterwards by
        :func:`organisms_released_after`.
    inputfile : str, optional
        Path to a file listing one tax id per line, used instead of ``taxon``
        when a download had to be expanded into many taxa to leave a branch out
        (see :mod:`eefinder.taxon_exclusion`). ``taxon`` is then only used for logging.

    Returns
    -------
    str
        The shell command (one token per space) to run.
    """
    if dataset not in DATASET_CHOICES:
        raise ValueError(f"Unknown dataset type: {dataset!r}")

    # A download restricted to a list of taxa names them in a file rather than
    # on the command line; the CLI takes one or the other, never both.
    selector = (
        f"--inputfile {shlex.quote(inputfile)}" if inputfile else shlex.quote(taxon)
    )
    if dataset == "virus":
        command = (
            f"{datasets_bin} download virus genome taxon {selector} "
            f"--include protein "
            f"--filename {zip_path}"
        )
        if refseq:
            command += " --refseq"
    else:  # bacteria / host: genome download, keep only the proteins
        command = (
            f"{datasets_bin} download genome taxon {selector} "
            f"--include protein "
            f"--filename {zip_path}"
        )
        if refseq:
            command += " --assembly-source RefSeq"
        if released_before:
            command += f" --released-before {shlex.quote(released_before)}"
    return command


class ConcatCounts(NamedTuple):
    """Result of :func:`concat_protein_faas`: file and sequence tallies."""

    files: int
    total: int
    written: int


def concat_protein_faas(
    extract_dir: str,
    output_fasta: str,
    exclude_products: "tuple[str, ...]" = (),
    tracking: "dict[str, dict] | None" = None,
    exclude_organisms: "set[str] | None" = None,
    release_dates: "dict[str, str] | None" = None,
) -> ConcatCounts:
    """Concatenate every ``protein.faa`` under ``extract_dir`` into one FASTA.

    A datasets download stores proteins under ``ncbi_dataset/data`` (one file
    for viruses, one per assembly for genomes), so all ``protein.faa`` files are
    merged in sorted order.

    Parameters
    ----------
    extract_dir : str
        Directory the datasets zip was extracted into.
    output_fasta : str
        Destination FASTA path.
    exclude_products : tuple[str, ...]
        Lower-cased product substrings whose records are dropped (e.g.
        :data:`UNINFORMATIVE_PRODUCTS`). Empty (the default) keeps everything.
    tracking : dict[str, dict], optional
        Filled with one row per **downloaded** accession, recording its product
        and whether it was dropped here. This is the only place the dropped ones
        are seen, since they never reach the FASTA.
    exclude_organisms : set[str], optional
        Organisms to leave out entirely (a release-date cutoff, see
        :func:`organisms_released_after`).
    release_dates : dict[str, str], optional
        ``organism -> earliest release date``, recorded for every accession so a
        dated build can be checked, not just trusted.

    Returns
    -------
    ConcatCounts
        ``files`` merged, ``total`` sequences seen, and ``written`` sequences
        kept (``total - written`` were dropped as uninformative).

    Raises
    ------
    FileNotFoundError
        If no ``protein.faa`` file is present in the download.
    """
    faas = sorted(Path(extract_dir).rglob("protein.faa"))
    if not faas:
        raise FileNotFoundError(f"no protein.faa found under {extract_dir}")
    exclude = tuple(term.lower() for term in exclude_products)
    total = written = 0
    keep = True
    with (
        open(output_fasta, "w") as out,
        progress_bar(faas, "Merging FASTAs ") as tracked,
    ):
        for faa in tracked:
            with open(faa) as handle:
                for line in handle:
                    if line.startswith(">"):
                        total += 1
                        header = parse_protein_header(line)
                        too_recent = bool(
                            exclude_organisms and header.organism in exclude_organisms
                        )
                        keep = not too_recent and not (
                            exclude
                            and any(term in header.product.lower() for term in exclude)
                        )
                        if keep:
                            written += 1
                        if tracking is not None:
                            tracking[header.accession] = {
                                "Accession": header.accession,
                                "Species": header.organism,
                                "Protein_downloaded": header.product,
                                "Protein_final": "",
                                "Name_changed": "",
                                "Status": STATUS_KEPT if keep else STATUS_REMOVED,
                                "Reason": (
                                    ""
                                    if keep
                                    else (
                                        REASON_TOO_RECENT
                                        if too_recent
                                        else REASON_UNINFORMATIVE
                                    )
                                ),
                                "Organism_release_date": (release_dates or {}).get(
                                    header.organism, ""
                                ),
                                "Cluster": "",
                                "Cluster_representative": "",
                            }
                    if keep:
                        out.write(line)
    return ConcatCounts(files=len(faas), total=total, written=written)


def count_fasta_records(fasta_path: str) -> int:
    """Return the number of records (header lines) in a FASTA file."""
    with open(fasta_path) as handle:
        return sum(1 for line in handle if line.startswith(">"))


def cluster_identical_proteins(
    fasta_path: str,
    threads: int = 1,
    clusters_path: "str | None" = None,
    tracking: "dict[str, dict] | None" = None,
) -> int:
    """Collapse 100%-identical / 100%-coverage duplicate proteins in place.

    Runs ``cd-hit`` at 100% sequence identity **and** 100% coverage of both the
    longer and the shorter sequence (``-c 1.0 -aL 1.0 -aS 1.0``), so only
    sequences that are identical over their entire length are merged to a single
    representative -- a lossless deduplication of the database. ``fasta_path`` is
    overwritten with the non-redundant set.

    Parameters
    ----------
    fasta_path : str
        Protein FASTA to deduplicate (overwritten in place).
    threads : int
        Threads for ``cd-hit`` (``-T``).
    clusters_path : str, optional
        Where to keep ``cd-hit``'s ``.clstr`` output, which lists the members of
        every cluster and marks the representative. Discarded when omitted.
    tracking : dict[str, dict], optional
        Updated with the cluster each accession landed in and the representative
        that stands for it, so a dropped duplicate can be traced to the sequence
        that replaced it.

    Returns
    -------
    int
        Number of duplicate sequences removed (``before - after``).
    """
    before = count_fasta_records(fasta_path)
    clustered = f"{fasta_path}.nr"
    cluster_proteins(fasta_path, clustered, threads)
    clstr = f"{clustered}.clstr"

    if tracking is not None and os.path.exists(clstr):
        for cluster_id, members in parse_cdhit_clusters(clstr).items():
            representative = members[0]
            for member in members:
                row = tracking.get(member)
                if row is None:
                    continue
                row["Cluster"] = cluster_id
                row["Cluster_representative"] = representative
                if member != representative:
                    row["Status"] = STATUS_REMOVED
                    row["Reason"] = REASON_DUPLICATE

    os.replace(clustered, fasta_path)
    if os.path.exists(clstr):
        if clusters_path:
            os.replace(clstr, clusters_path)
        else:
            os.remove(clstr)
    after = count_fasta_records(fasta_path)
    return before - after


def validate_date(value: "str | None") -> "str | None":
    """Check a ``YYYY-MM-DD`` cutoff, raising a clear error when it is not one.

    Returns ``None`` unchanged, so callers can pass an unset option straight
    through.
    """
    if value is None or value == "":
        return None
    try:
        datetime.strptime(value, "%Y-%m-%d")
    except ValueError:
        raise ValueError(f"Invalid date {value!r}: use ISO format, e.g. 2024-12-31.")
    return value


def reconcile_tracking(tracking: "dict[str, dict]", fasta_path: str) -> int:
    """Mark as removed anything still called ``kept`` that is not in the FASTA.

    The tracking table must agree with the database it describes. Steps can drop
    a record without saying so -- ``cd-hit`` silently discards sequences of at
    most :data:`CDHIT_MIN_LENGTH` residues -- and such a record would otherwise
    be reported as kept while being absent from both the FASTA and the CSV.

    Parameters
    ----------
    tracking : dict[str, dict]
        The table built through the run.
    fasta_path : str
        The final protein FASTA.

    Returns
    -------
    int
        Number of rows corrected.
    """
    present = set()
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                present.add(line[1:].split(None, 1)[0])

    corrected = 0
    for accession, row in tracking.items():
        if row["Status"] == STATUS_KEPT and accession not in present:
            row["Status"] = STATUS_REMOVED
            row["Reason"] = REASON_ABSENT
            corrected += 1
    return corrected


def write_tracking_table(tracking: "dict[str, dict]", out_path: str) -> int:
    """Write the per-accession audit table (:data:`TRACKING_COLUMNS`).

    One row per sequence **as downloaded**, so every accession that did not make
    it into the final database can be traced to the reason it left: an
    uninformative product, an identical duplicate (with the cluster and the
    representative that replaced it), or a product that standardised to
    ``Unknown``. Kept sequences record whether standardisation renamed them.

    Returns
    -------
    int
        Number of rows written.
    """
    frame = pd.DataFrame(list(tracking.values()), columns=TRACKING_COLUMNS)
    frame = frame.sort_values(by=["Status", "Accession"])
    frame.to_csv(out_path, sep="\t", index=False)
    return len(frame)


def find_data_reports(extract_dir: str) -> "list[str]":
    """Return every ``data_report.jsonl`` under ``extract_dir``.

    A download split into several packages (see :mod:`eefinder.taxon_exclusion`) is
    extracted one package per subdirectory, so there is a report per part.
    """
    return [str(path) for path in sorted(Path(extract_dir).rglob("data_report.jsonl"))]


def find_data_report(extract_dir: str) -> str | None:
    """Return the first ``data_report.jsonl`` under ``extract_dir``, if any."""
    reports = find_data_reports(extract_dir)
    return reports[0] if reports else None


def filter_fasta_by_ids(fasta_path: str, keep_ids: "set[str]") -> int:
    """Rewrite ``fasta_path`` in place, keeping only records in ``keep_ids``.

    A record is kept when the first whitespace-delimited token of its header
    (i.e. its accession/id) is in ``keep_ids``. Used to drop from the FASTA the
    same records dropped from the metadata CSV, keeping the two in sync.

    Parameters
    ----------
    fasta_path : str
        FASTA to filter (overwritten via a temporary file).
    keep_ids : set[str]
        Header ids (first token, without ``>``) to retain.

    Returns
    -------
    int
        Number of records dropped.
    """
    tmp_path = f"{fasta_path}.tmp"
    dropped = 0
    keep = True
    with open(fasta_path) as fasta_in, open(tmp_path, "w") as fasta_out:
        for line in fasta_in:
            if line.startswith(">"):
                record_id = line[1:].split(None, 1)[0]
                keep = record_id in keep_ids
                if not keep:
                    dropped += 1
            if keep:
                fasta_out.write(line)
    os.replace(tmp_path, fasta_path)
    return dropped


class GetDatabases:
    """Download and assemble an EEfinder database via the NCBI datasets CLI.

    Runs on instantiation: downloads the datasets zip, extracts it, writes the
    concatenated protein FASTA (``{outdir}/{prefix}.fa``) and, for the ``virus``
    and ``bacteria`` datasets, the metadata CSV (``{outdir}/{prefix}.csv``). A
    JSON run summary with kept/dropped sequence counts is written to
    ``{outdir}/{prefix}.log`` and exposed as ``self.sequence_counts``.

    Parameters
    ----------
    dataset : str
        One of :data:`DATASET_CHOICES`.
    taxon : str
        NCBI taxon name or tax id to download.
    outdir : str
        Output directory (created if missing).
    prefix : str, optional
        Basename for the output files; defaults to ``dataset``.
    refseq : bool
        Restrict to RefSeq sequences (recommended).
    exclude_uninformative : bool
        Drop ``hypothetical protein`` / ``uncharacterized protein`` records
        (:data:`UNINFORMATIVE_PRODUCTS`) from the FASTA and CSV.
    standardize_proteins : bool
        Rewrite the CSV ``Protein`` column to canonical names via the bundled
        viral protein map (:func:`standardize_protein`).
    cluster : bool
        Collapse 100%-identical / 100%-coverage duplicate proteins with
        ``cd-hit`` (:func:`cluster_identical_proteins`) before writing the
        metadata CSV.
    threads : int
        Threads for the ``cd-hit`` clustering step.
    datasets_bin : str
        Path/name of the ``datasets`` executable.
    attempts : int
        How many times the download is attempted before giving up.
    stall_timeout : float
        Seconds without any sign of life (no output from ``datasets``, no growth
        of the archive) after which an attempt is treated as hung, killed and
        retried.
    keep_download : bool
        Keep the downloaded zip and the directory it was extracted into. They
        are deleted by default: their content is already in the outputs.
    released_before : str, optional
        Only include data released on or before this date (``YYYY-MM-DD``), so a
        build can be reproduced later. Applied by the ``datasets`` CLI for
        bacteria/host and by EEfinder for viruses, whose subcommand has no such
        flag.
    exclude_taxa : tuple[str, ...]
        Tax ids or scientific names to leave out of the download. The branch is
        never requested: ``taxon`` is expanded into the taxa that cover it minus
        the excluded subtrees (:func:`eefinder.taxon_exclusion.expand_taxon_excluding`),
        which is what keeps SARS-CoV-2 -- 61% of all viral records -- out of an
        ``--all-sequences`` viral build.
    """

    def __init__(
        self,
        dataset: str,
        taxon: str,
        outdir: str,
        prefix: "str | None" = None,
        refseq: bool = True,
        exclude_uninformative: bool = False,
        standardize_proteins: bool = False,
        cluster: bool = False,
        threads: int = 1,
        datasets_bin: str = DATASETS_BINARY,
        attempts: int = DEFAULT_ATTEMPTS,
        stall_timeout: float = DEFAULT_STALL_TIMEOUT,
        keep_download: bool = False,
        released_before: "str | None" = None,
        exclude_taxa: "tuple[str, ...]" = (),
    ) -> None:
        if dataset not in DATASET_CHOICES:
            raise ValueError(f"Unknown dataset type: {dataset!r}")
        self.dataset = dataset
        self.taxon = taxon
        self.outdir = check_outdir(outdir)
        self.prefix = prefix or dataset
        self.refseq = refseq
        self.exclude_uninformative = exclude_uninformative
        self.standardize_proteins = standardize_proteins
        self.cluster = cluster
        self.threads = threads
        self.datasets_bin = datasets_bin
        self.attempts = attempts
        self.stall_timeout = stall_timeout
        self.keep_download = keep_download
        self.released_before = validate_date(released_before)
        self.exclude_taxa = tuple(exclude_taxa or ())

        self.get_databases()

    def get_databases(self) -> None:
        """Download the dataset and write the FASTA (and CSV where applicable)."""
        run_start = time.time()
        steps: list[StepInfo] = []
        tracking: "dict[str, dict]" = {}
        zip_path = f"{self.outdir}/{self.prefix}.zip"
        extract_dir = f"{self.outdir}/{self.prefix}_ncbi"
        fasta_out = f"{self.outdir}/{self.prefix}.fa"
        logger.debug(
            f"Paths: zip={zip_path} extract_dir={extract_dir} fasta={fasta_out}"
        )

        start = time.time()
        batches = self._plan_download()
        if len(batches) > 1 or (batches and batches[0] != [self.taxon]):
            steps.append(
                StepInfo.from_times(
                    "Resolve taxonomy",
                    start,
                    time.time(),
                    self._exclusion_summary,
                )
            )

        logger.info(f"Downloading {self.dataset} proteins for taxon '{self.taxon}'")
        start = time.time()
        zip_paths = self._download_batches(batches, zip_path)
        steps.append(
            StepInfo.from_times(
                "Download",
                start,
                time.time(),
                f"Downloaded {self.dataset} proteins for taxon '{self.taxon}' "
                f"(refseq={self.refseq}) to "
                f"{', '.join(zip_paths) if len(zip_paths) > 1 else zip_path}.",
            )
        )

        logger.info("Extracting the datasets archive")
        start = time.time()
        for index, part_zip in enumerate(zip_paths, start=1):
            # One package per subdirectory: every package carries its own
            # data_report.jsonl and protein.faa, which would otherwise overwrite
            # each other.
            target = (
                extract_dir if len(zip_paths) == 1 else f"{extract_dir}/part_{index}"
            )
            with zipfile.ZipFile(part_zip) as archive:
                members = archive.infolist()
                label = (
                    "Extracting     "
                    if len(zip_paths) == 1
                    else f"Extracting {index}/{len(zip_paths)}"
                )
                with progress_bar(members, label) as tracked:
                    for member in tracked:
                        archive.extract(member, target)
        logger.debug(f"Extracted archive into {extract_dir}")
        steps.append(
            StepInfo.from_times(
                "Extract",
                start,
                time.time(),
                f"Extracted {len(zip_paths)} package(s) into {extract_dir}.",
            )
        )

        exclude = UNINFORMATIVE_PRODUCTS if self.exclude_uninformative else ()
        too_recent: "set[str]" = set()
        release_dates: "dict[str, str]" = {}
        reports = find_data_reports(extract_dir)
        for report in reports:
            # Recorded for every accession, whether or not a cutoff was asked
            # for: a dated build should be checkable, not merely trusted.
            release_dates.update(parse_release_dates(report))
        if self.released_before and self.dataset == "virus":
            # The virus subcommand of datasets has no --released-before, so the
            # cutoff is applied here, from the release dates in the report.
            if reports:
                too_recent = {
                    organism
                    for organism, released in release_dates.items()
                    if released > self.released_before
                }
                logger.info(
                    f"Release cutoff {self.released_before}: leaving out "
                    f"{len(too_recent)} organism(s) that first appeared later"
                )
            else:
                logger.warning(
                    "No data_report.jsonl in the download; the release cutoff "
                    "could not be applied."
                )
        start = time.time()
        counts = concat_protein_faas(
            extract_dir,
            fasta_out,
            exclude_products=exclude,
            tracking=tracking,
            exclude_organisms=too_recent,
            release_dates=release_dates,
        )
        excluded_uninformative = counts.total - counts.written
        merged = f"{counts.files} protein.faa file(s) merged"
        if exclude:
            merged += "; dropped hypothetical/uncharacterized proteins"
        logger.info(f"Wrote {fasta_out} ({merged})")
        logger.debug(
            f"Sequences: downloaded={counts.total} "
            f"excluded_uninformative={excluded_uninformative} "
            f"written={counts.written}"
        )
        steps.append(
            StepInfo.from_times(
                "Write protein FASTA",
                start,
                time.time(),
                f"Merged {counts.files} protein.faa file(s): {counts.total} "
                f"sequences, {excluded_uninformative} dropped as uninformative, "
                f"{counts.written} written to {fasta_out}.",
            )
        )

        clustered_identical = 0
        if self.cluster:
            start = time.time()
            clustered_identical = cluster_identical_proteins(
                fasta_out,
                self.threads,
                clusters_path=f"{self.outdir}/{self.prefix}.clstr",
                tracking=tracking,
            )
            remaining = counts.written - clustered_identical
            logger.info(
                f"Clustered {fasta_out}: removed {clustered_identical} "
                f"100%/100% duplicate(s), {remaining} representative(s) kept"
            )
            steps.append(
                StepInfo.from_times(
                    "Cluster identical proteins",
                    start,
                    time.time(),
                    f"cd-hit (100% identity / 100% coverage): removed "
                    f"{clustered_identical} duplicate(s), {remaining} of "
                    f"{counts.written} representative(s) kept in {fasta_out}.",
                )
            )

        dropped_standardization = 0
        if self.dataset in _METADATA_DATASETS:
            start = time.time()
            reports = find_data_reports(extract_dir)
            logger.debug(f"data_report.jsonl: {reports}")
            taxonomy: "dict[str, TaxonomyRecord]" = {}
            for report in reports:
                taxonomy.update(parse_taxonomy_report(report))
            logger.debug(f"Parsed taxonomy for {len(taxonomy)} organism(s)")
            if not reports:
                logger.warning(
                    "No data_report.jsonl in the download; the metadata CSV will "
                    "have empty Genus/Family/Molecule_type/Host columns."
                )
            frame = build_metadata_frame(
                fasta_out,
                taxonomy,
                standardize=self.standardize_proteins,
                dataset=self.dataset,
                tracking=tracking,
            )
            if self.standardize_proteins:
                # Keep the FASTA in sync: drop the records that standardisation
                # removed from the table (bare CDS/ORF and hypothetical proteins).
                dropped_standardization = filter_fasta_by_ids(
                    fasta_out, set(frame["Accession"])
                )
                if dropped_standardization:
                    logger.info(
                        f"Dropped {dropped_standardization} unknown "
                        f"(CDS/ORF/hypothetical) protein(s) from {fasta_out}"
                    )
            csv_out = f"{self.outdir}/{self.prefix}.csv"
            frame.to_csv(csv_out, index=False)
            logger.info(f"Wrote {csv_out} ({len(frame)} records)")
            steps.append(
                StepInfo.from_times(
                    "Write metadata CSV",
                    start,
                    time.time(),
                    f"Wrote {len(frame)} records to {csv_out} "
                    f"(standardize={self.standardize_proteins}); dropped "
                    f"{dropped_standardization} sequence(s) from the FASTA.",
                )
            )

        start = time.time()
        unexplained = reconcile_tracking(tracking, fasta_out)
        if unexplained:
            logger.warning(
                f"{unexplained} sequence(s) are absent from {fasta_out} without a "
                "step reporting it; they are marked "
                f"'{REASON_ABSENT}' in the tracking table"
                + (
                    f" (cd-hit discards sequences of up to {CDHIT_MIN_LENGTH} "
                    "residues)"
                    if self.cluster
                    else ""
                )
            )
        tracking_path = f"{self.outdir}/{self.prefix}.tracking.tsv"
        rows = write_tracking_table(tracking, tracking_path)
        logger.info(f"Wrote {tracking_path} ({rows} downloaded accession(s))")
        steps.append(
            StepInfo.from_times(
                "Write tracking table",
                start,
                time.time(),
                f"Recorded the fate of {rows} downloaded accession(s) in "
                f"{tracking_path}.",
            )
        )

        if not self.keep_download:
            start = time.time()
            removed = self._clean_download(zip_path, extract_dir)
            steps.append(
                StepInfo.from_times(
                    "Remove the raw download",
                    start,
                    time.time(),
                    f"Deleted {removed} (the database, metadata, tracking table "
                    "and logs are kept).",
                )
            )

        kept = counts.written - clustered_identical - dropped_standardization
        self.sequence_counts = SequenceCounts(
            downloaded=counts.total,
            excluded_uninformative=excluded_uninformative,
            clustered_identical=clustered_identical,
            dropped_standardization=dropped_standardization,
            kept=kept,
        )
        logger.info(
            f"Sequences: {kept} kept, "
            f"{excluded_uninformative + clustered_identical + dropped_standardization}"
            f" dropped (of {counts.total} downloaded)"
        )
        self._write_log(run_start, time.time(), steps)

    def _clean_download(self, zip_path: str, extract_dir: str) -> str:
        """Delete the downloaded archive and the directory it was extracted to.

        Everything they contain has been read into the outputs by this point;
        for a whole-RefSeq download they are several GB of duplication. What is
        kept is the database FASTA, the metadata CSV, the tracking table and the
        logs.
        """
        removed = []
        base = zip_path[: -len(".zip")] if zip_path.endswith(".zip") else zip_path
        # A download split across taxa leaves {prefix}.partN.zip beside the
        # single-package {prefix}.zip; both spellings are cleaned up.
        archives = [zip_path] + sorted(glob.glob(f"{base}.part*.zip"))
        for archive in archives:
            if os.path.exists(archive):
                os.remove(archive)
                removed.append(os.path.basename(archive))
        if os.path.isdir(extract_dir):
            shutil.rmtree(extract_dir, ignore_errors=True)
            removed.append(f"{os.path.basename(extract_dir)}/")
        if removed:
            logger.info(f"Removed the raw download: {', '.join(removed)}")
        return ", ".join(removed) or "nothing"

    def _write_log(
        self, start_time: float, end_time: float, steps: list[StepInfo]
    ) -> None:
        """Write the ``{prefix}.log`` JSON run summary (kept/dropped counts)."""
        info = DownloadInfo.from_run(
            eefinder_version=__version__,
            arguments=DownloadArguments(
                dataset=self.dataset,
                taxon=self.taxon,
                outdir=self.outdir,
                prefix=self.prefix,
                refseq=self.refseq,
                exclude_uninformative=self.exclude_uninformative,
                standardize_proteins=self.standardize_proteins,
                cluster=self.cluster,
                attempts=self.attempts,
                stall_timeout=self.stall_timeout,
                keep_download=self.keep_download,
                released_before=self.released_before or "",
                exclude_taxa=", ".join(self.exclude_taxa),
            ),
            sequence_counts=self.sequence_counts,
            start_time=start_time,
            end_time=end_time,
            steps_information=steps,
        )
        log_path = f"{self.outdir}/{self.prefix}.log"
        logger.debug(f"Writing download summary to {log_path}")
        with open(log_path, "w") as json_out:
            json.dump(asdict(info), json_out, indent=4)
        logger.info(f"Wrote {log_path}")

    def _plan_download(self) -> "list[list[str]]":
        """Work out which taxa to ask for, honouring ``exclude_taxa``.

        Returns one list of taxa per ``datasets`` call: ``[[self.taxon]]`` when
        nothing is excluded, otherwise the sibling taxa that cover the requested
        taxon minus the excluded branches, split into ``--inputfile``-sized
        batches. The excluded branch is never requested, so its records are
        never transferred.
        """
        self._exclusion_summary = ""
        if not self.exclude_taxa:
            return [[self.taxon]]

        logger.info("Resolving taxonomy to leave out: " + ", ".join(self.exclude_taxa))
        expansion = expand_taxon_excluding(
            self.taxon, self.exclude_taxa, datasets_bin=self.datasets_bin
        )
        for node in expansion.not_below_root:
            logger.warning(
                f"{node.name} ({node.tax_id}) is not below "
                f"{expansion.root.name} ({expansion.root.tax_id}); nothing to "
                "exclude for it"
            )
        if not expansion.pruned:
            return [[self.taxon]]

        left_out = ", ".join(
            f"{node.name} ({node.tax_id})" for node in expansion.excluded
        )
        batches = batch_taxa(expansion.taxa)
        self._exclusion_summary = (
            f"Expanded taxon '{self.taxon}' into {len(expansion.taxa)} taxa in "
            f"{len(batches)} download(s) to leave out {left_out}. Records "
            "attached directly to a rank between them are not reachable this way."
        )
        logger.info(
            f"Leaving out {left_out}: downloading {len(expansion.taxa)} sibling "
            f"taxa instead of '{self.taxon}'"
        )
        return [[str(t) for t in batch] for batch in batches]

    def _download_batches(
        self, batches: "list[list[str]]", zip_path: str
    ) -> "list[str]":
        """Download each batch of taxa, returning the packages written."""
        if len(batches) == 1 and batches[0] == [self.taxon]:
            self._download(zip_path)
            return [zip_path]

        base = zip_path[: -len(".zip")] if zip_path.endswith(".zip") else zip_path
        written: "list[str]" = []
        for index, batch in enumerate(batches, start=1):
            part_zip = zip_path if len(batches) == 1 else f"{base}.part{index}.zip"
            listing = f"{base}.part{index}.taxa.txt"
            with open(listing, "w") as handle:
                handle.write("\n".join(batch) + "\n")
            logger.info(
                f"Downloading package {index}/{len(batches)} ({len(batch)} taxa)"
            )
            try:
                self._download(part_zip, inputfile=listing)
            finally:
                if not self.keep_download and os.path.exists(listing):
                    os.remove(listing)
            written.append(part_zip)
        return written

    def _download(self, zip_path: str, inputfile: "str | None" = None) -> None:
        """Run ``datasets download``, raising on a non-zero exit."""
        command = build_download_command(
            self.dataset,
            self.taxon,
            zip_path,
            refseq=self.refseq,
            datasets_bin=self.datasets_bin,
            inputfile=inputfile,
        )
        logger.debug(f"datasets command: {command}")

        def _archive_is_complete() -> bool:
            """Whether the download already produced a readable, complete zip.

            ``zipfile`` can only read the central directory once the archive has
            been written in full, so this cannot be satisfied by a partial file.
            It is what lets the run continue when ``datasets`` has delivered the
            package but does not exit.
            """
            if not os.path.exists(zip_path) or not zipfile.is_zipfile(zip_path):
                return False
            try:
                with zipfile.ZipFile(zip_path) as archive:
                    return bool(archive.namelist()) and archive.testzip() is None
            except (zipfile.BadZipFile, OSError):
                return False

        def _discard_partial(attempt: "int | None" = None) -> None:
            """Drop any archive already at ``zip_path`` so a run starts clean.

            Called before the first attempt as well as between retries: a
            leftover archive from an interrupted run is not something to download
            *into*.
            """
            if os.path.exists(zip_path):
                logger.debug(f"removing pre-existing download {zip_path}")
                os.remove(zip_path)

        _discard_partial()

        def _http1_fallback(output: str) -> "dict[str, str] | None":
            """Retry over HTTP/1.1 after NCBI resets an HTTP/2 stream.

            ``stream error: stream ID 3; INTERNAL_ERROR; received from peer``
            aborts large transfers part-way through; the same download completes
            over HTTP/1.1, so repeating the attempt unchanged would just fail the
            same way.
            """
            lowered = (output or "").lower()
            if any(marker in lowered for marker in HTTP2_ERROR_MARKERS):
                logger.warning(
                    "NCBI reset the HTTP/2 stream; the retry will use HTTP/1.1"
                )
                return dict(HTTP1_ENV)
            return None

        # The datasets CLI draws its own download/validation progress; mirror it
        # to the terminal instead of swallowing it, while keeping the tail so a
        # failure can still be reported with the tool's own message. NCBI
        # transfers also hang outright, which no exit status reports, so the run
        # is watched for a stall and retried.
        returncode, output = run_with_retries(
            shlex.split(command),
            attempts=self.attempts,
            stall_timeout=self.stall_timeout,
            liveness_path=zip_path,
            before_retry=_discard_partial,
            success_check=_archive_is_complete,
            retry_env=_http1_fallback,
        )
        if returncode == 0 and not _archive_is_complete():
            # datasets can report success and still leave a broken package.
            returncode, output = STALLED, (
                output or "the download finished but the archive is not readable"
            )
            logger.warning(
                f"{zip_path} is not a readable zip archive despite a clean exit"
            )

        if returncode != 0:
            message = output
            if "ACELLULAR_ROOT" in message or "RankType" in message:
                message += (
                    " -- this is an outdated NCBI datasets CLI: it predates the "
                    "'acellular root' taxonomy rank added above Viruses. Upgrade "
                    "to datasets >= 18.1 (see env.yml)."
                )
            reason = (
                f"it stalled for more than {self.stall_timeout:.0f}s"
                if returncode == STALLED
                else f"exit {returncode}"
            )
            raise RuntimeError(
                f"datasets download failed after {self.attempts} attempt(s) "
                f"({reason}): {message}"
            )
