#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Command-line entrypoint orchestrating the EEfinder pipeline.

Each pipeline step is a small class in :mod:`eefinder` that runs its work on
instantiation and communicates with the next step through files on disk. This
module wires those steps together in order and writes a typed run summary
(:class:`~eefinder.utils.RunInfo`) to ``eefinder.log``.
"""

import faulthandler
import io
import signal
import time

import click
import re
import os
import glob
import shutil
import sys
import json
from dataclasses import asdict
from eefinder.log import logger, enable_debug
from eefinder.run_message import PaperInfo
from eefinder.utils import check_outdir, StepInfo, RunArguments, RunInfo
from eefinder.prepare_data import PrepareGenome
from eefinder.clean_data import MaskClean
from eefinder.make_database import MakeDB
from eefinder.similarity_analysis import SimilaritySearch
from eefinder.translation import TRANSLATION_METHODS
from eefinder.filter_table import FilterTable
from eefinder.get_taxonomy import GetTaxonomy, GetFinalTaxonomy, GetCleanedTaxonomy
from eefinder.bed import (
    GetFasta,
    GetAnnotBed,
    RemoveAnnotation,
    MergeBed,
    BedFlank,
    GetBed,
)
from eefinder.compare_results import CompareResults
from eefinder.get_length import GetLength
from eefinder.tag_elements import TagElements
from eefinder.overlap import FilterOverlap
from eefinder.get_databases import (
    GetDatabases,
    DATASETS_BINARY,
    CDHIT_BINARY,
    DEFAULT_TAXA,
)
from eefinder.taxon_exclusion import DEFAULT_VIRUS_EXCLUSIONS, NO_EXCLUSION
from eefinder.progress import DEFAULT_ATTEMPTS, DEFAULT_STALL_TIMEOUT
from eefinder.gff import WriteGFF3
from eefinder.versions import (
    collect_dependency_versions,
    collect_system_info,
    find_env_yml,
)
from eefinder import __version__

HOMEPAGE = "https://github.com/WallauBioinfo/EEfinder"


def _enable_stack_dumps() -> None:
    """Let ``kill -USR1 <pid>`` print what every thread is doing.

    A run that appears frozen is otherwise impossible to diagnose after the
    fact: this turns it into one command that says exactly which line is
    waiting, including inside external-tool handling.
    """
    if not hasattr(signal, "SIGUSR1"):  # not available on Windows
        return
    try:
        # sys.__stderr__ rather than sys.stderr: the handler needs a real file
        # descriptor, which an in-memory replacement (a test runner's, say)
        # does not have.
        faulthandler.register(
            signal.SIGUSR1, file=sys.__stderr__, all_threads=True, chain=True
        )
    except (AttributeError, ValueError, OSError, io.UnsupportedOperation):
        logger.debug("stack dumps on SIGUSR1 are not available here")


def report_run_context(system, dependencies):
    """Log the EEfinder version, host context and dependency versions.

    Prokka-style startup banner. Emits a warning for any dependency whose
    runtime version differs from the ``env.yml`` pin or that could not be found
    on ``PATH``.
    """
    logger.info(f"This is EEfinder {__version__}")
    logger.info(f"Homepage is {HOMEPAGE}")
    logger.info(f"Operating system is {system.operating_system}")
    logger.info(f"You are {system.user}")
    for dep in dependencies:
        if dep.status == "mismatch":
            logger.warning(
                f"{dep.name} {dep.detected} differs from the env.yml pin "
                f"{dep.expected}"
            )
        elif dep.status == "not-found":
            logger.warning(f"{dep.name} was not found on PATH")
        else:
            logger.info(f"Using {dep.name} {dep.detected}")


@click.group()
@click.version_option(__version__)
def cli():
    """EEfinder: find Endogenous Elements in eukaryote genomes.

    Use ``screening`` to run the EE-finding pipeline and ``get-databases`` to
    download the RefSeq protein databases it needs.
    """
    _enable_stack_dumps()


@cli.command(name="screening")
@click.version_option(__version__)
@click.option(
    "-in",
    "--genome_file",
    help="Input genome fasta file (nucleotides).",
    required=True,
)
@click.option(
    "-od",
    "--outdir",
    help="Path and dir to store output results.",
    required=True,
)
@click.option(
    "-db",
    "--database",
    help="Proteins from viruses or bacterias database .fasta file.",
    required=True,
)
@click.option(
    "-mt",
    "--dbmetadata",
    help="Proteins from viruses or bacterias metadata .csv file.",
    required=True,
)
@click.option(
    "-bt",
    "--hostgenesbaits",
    help="Host genes baits proteins, used to filter putative EEs .fasta file.",
    required=True,
)
@click.option(
    "-md",
    "--mode",
    help="Choose between BLAST or the DIAMOND strategies (fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensisitve) to run analysis, default = blastx.",
    default="blastx",
    type=click.Choice(
        [
            "blastx",
            "fast",
            "mid-sensitive",
            "sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ]
    ),
)
@click.option(
    "-ln",
    "--length",
    help="Minimum length of contigs used for BLAST or DIAMOND, default = 10000.",
    type=int,
    default=10000,
)
@click.option(
    "-fl",
    "--flank",
    help="Length of flanking regions of Endogenous Elements to be extracted, default = 10000.",
    type=int,
    default=10000,
)
@click.option(
    "-lm",
    "--limit",
    help="Limit of bases used to merge regions on bedtools merge, default = 1.",
    type=int,
    default=1,
)
@click.option(
    "-rj",
    "--range_junction",
    help="Sets the range for junction of BLAST/DIAMOND redudant hits, default=100",
    type=int,
    default=100,
)
@click.option(
    "-mp",
    "--mask_per",
    help="Limit of lowercase letters in percentage to consider a putative Endogenous Elements as a repetitive region, default = 50.",
    type=int,
    default=50,
)
@click.option(
    "-cm",
    "--clean_masked",
    help="Remove EEs in regions considered repetitive?",
    is_flag=True,
)
@click.option(
    "-p",
    "--threads",
    help="Threads for multi-thread analysis, default = 1.",
    type=int,
    default=1,
)
@click.option(
    "-rm",
    "--removetmp",
    help="Remove temporary files generated through analysis?",
    is_flag=True,
)
@click.option(
    "-id",
    "--index_databases",
    help="Index databases?",
    is_flag=True,
)
@click.option(
    "-pr",
    "--prefix",
    help="Write the prefix name for output files. This prefix will be used to create the EEname (The Endogenous Element name will be formated as PREFIX|CONTIG/SCAFFOLD:START-END) default = input file name.",
)
@click.option(
    "-ml",
    "--merge_level",
    help="Taxonomy level to merge elements by genus or family, default = family",
    default="family",
    type=click.Choice(["family", "genus"]),
)
@click.option(
    "-an",
    "--analysis",
    help="Type of endogenous elements being screened; sets the GFF3 feature "
    "type ('virus' -> endogenous_viral_element, 'bacteria' -> "
    "endogenous_bacterial_element). default = virus",
    default="virus",
    type=click.Choice(["virus", "bacteria"]),
)
@click.option(
    "-ov",
    "--overlap",
    help="What to do with elements tagged as overlaped: 'keep' (default, keep "
    "all), 'longest' (keep the longest element of each overlap), or 'targets' "
    "(keep only elements from --target_families). Filtered-out elements are "
    "saved to tmp_outputs/. default = keep",
    default="keep",
    type=click.Choice(["keep", "longest", "targets"]),
)
@click.option(
    "-tf",
    "--target_families",
    help="Family to KEEP when --overlap targets is used. Repeat for multiple "
    "families (e.g. -tf Flaviviridae -tf Caulimoviridae). Mutually exclusive "
    "with --non_target_families.",
    multiple=True,
)
@click.option(
    "-ntf",
    "--non_target_families",
    help="Family to DROP when --overlap targets is used. Repeat for multiple "
    "families (e.g. -ntf Retroviridae -ntf Metaviridae). Mutually exclusive "
    "with --target_families.",
    multiple=True,
)
@click.option(
    "-tm",
    "--translation_method",
    help="How proteins are obtained for the similarity searches (applied to "
    "BOTH the main and the host-bait search): 'default' = six-frame "
    "blastx/diamond blastx; 'gv' = pyrodigal-gv prediction; 'rv' = pyrodigal-rv "
    "prediction; 'gv-rv' = both predictions clustered with cd-hit (100%/100%). "
    "Prediction modes align with blastp and map coordinates back to nucleotides. "
    "default = default",
    default="default",
    type=click.Choice(list(TRANSLATION_METHODS)),
)
@click.option(
    "--debug",
    help="Emit verbose debug logging (intermediate file paths, per-step "
    "details). default = off",
    is_flag=True,
)
def screening(
    genome_file,
    outdir,
    database,
    dbmetadata,
    hostgenesbaits,
    mode,
    length,
    flank,
    limit,
    range_junction,
    mask_per,
    clean_masked,
    threads,
    removetmp,
    index_databases,
    prefix,
    merge_level,
    analysis,
    overlap,
    target_families,
    non_target_families,
    translation_method,
    debug,
):
    """Run the EEfinder screening pipeline on a genome."""
    if debug:
        enable_debug()
    steps_infos: "list[StepInfo]" = []
    start_running_time = time.time()
    print_info = PaperInfo()
    print_info.print_start(__version__)
    logger.debug(
        "screening arguments: "
        f"genome_file={genome_file!r} outdir={outdir!r} database={database!r} "
        f"dbmetadata={dbmetadata!r} baits={hostgenesbaits!r} mode={mode!r} "
        f"length={length} flank={flank} limit={limit} "
        f"range_junction={range_junction} mask_per={mask_per} "
        f"clean_masked={clean_masked} threads={threads} "
        f"index_databases={index_databases} prefix={prefix!r} "
        f"merge_level={merge_level!r} analysis={analysis!r} overlap={overlap!r} "
        f"target_families={list(target_families)} "
        f"non_target_families={list(non_target_families)} "
        f"translation_method={translation_method!r}"
    )

    system_info = collect_system_info()
    try:
        dependencies = collect_dependency_versions(find_env_yml())
        report_run_context(system_info, dependencies)
    except Exception as err:
        logger.warning(f"Could not collect dependency versions: {err}")
        dependencies = []

    if overlap == "targets" and bool(target_families) == bool(non_target_families):
        click.secho(
            "--overlap targets requires exactly one of --target_families or "
            "--non_target_families.",
            err=True,
            fg="red",
        )
        sys.exit(1)

    if prefix is None:
        try:
            logger.info(f"Creating prefix")
            prefix = genome_file
            prefix = re.sub(r"\..*", "", prefix)
            prefix = re.sub(r".*/", "", prefix).rstrip("\n")
            prefix = re.sub(r".*/", "", prefix).rstrip("\n")
        except Exception as err:
            click.secho(f"Failed to create prefix: {err}", err=True, fg="red")
            sys.exit(1)
    logger.debug(f"Using prefix {prefix!r}")

    try:
        logger.info(f"Creating output directory")
        outdir = check_outdir(outdir)
        logger.debug(f"Output directory resolved to {outdir}")

    except Exception as err:
        click.secho(f"Failed to create output dir: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info(f"Preparing input data")
        start_time = time.time()

        logger.debug(
            f"PrepareGenome: {genome_file} -> {outdir}/{prefix}.rn.fmt "
            f"(prefixing headers, dropping contigs < {length} nt)"
        )
        prepared = PrepareGenome(genome_file, prefix, outdir, length)
        logger.info(
            f"Prepared {prepared.kept} of {prepared.total} contig(s) "
            f"(>= {length} nt)"
        )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Prepare input data",
                start_time=start_time,
                end_time=end_time,
                message=f"{prefix} prefix included in {genome_file} sequences header and sequences bellow than {length} nt are removed from {genome_file}.",
            )
        )
    except Exception as err:
        click.secho(f"Failed to prepare input data: {err}", err=True, fg="red")
        sys.exit(1)

    if index_databases:
        try:
            logger.info(f"Indexing databases")
            start_time = time.time()

            logger.debug(f"MakeDB ({mode}) protein DB: {database}")
            MakeDB(mode, database, "prot", threads)
            logger.debug(f"MakeDB ({mode}) baits DB: {hostgenesbaits}")
            MakeDB(mode, hostgenesbaits, "prot", threads)

            end_time = time.time()
            steps_infos.append(
                StepInfo.from_times(
                    step="Index databases",
                    start_time=start_time,
                    end_time=end_time,
                    message=f"Index {database} and {hostgenesbaits} databases.",
                )
            )
        except Exception as err:
            click.secho(f"Failed to format databases: {err}", err=True, fg="red")
            sys.exit(1)
    else:
        logger.warning("index_databases step will not be performed")

    try:
        logger.info(f"Running similarity search")
        start_time = time.time()
        query = f"{outdir}/{prefix}.rn.fmt"

        logger.debug(
            f"SimilaritySearch ({mode}, translation={translation_method}): "
            f"{query} vs {database} -> {query}.blastx"
        )
        SimilaritySearch(query, database, threads, mode, translation_method)
        logger.debug(
            f"FilterTable EE: {query}.blastx (range_junction={range_junction})"
        )
        FilterTable(f"{query}.blastx", range_junction, "EE", outdir)

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Similarity search",
                start_time=start_time,
                end_time=end_time,
                message=f"Similarity analysis with {mode} was performed using {query} against {database}."
                + f"Matches against same subject sequence in a {range_junction}nt range junction are filtered, mantaining the one with the greatest bitscore.",
            )
        )
    except Exception as err:
        click.secho(f"Failed to run similarity searches: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Extracting Putative EEs")
        start_time = time.time()

        logger.debug(
            f"GetFasta putative EEs from {outdir}/{prefix}.rn.fmt.blastx.filtred.bed"
        )
        GetFasta(
            f"{outdir}/{prefix}.rn.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta",
        )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Extraction of putative EEs",
                start_time=start_time,
                end_time=end_time,
                message=f"Extract infos based on {outdir}/{prefix}.rn.fmt.blastx.filtred.bed information.",
            )
        )
    except Exception as err:
        click.secho(f"Failed to run extract EE sequences: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Running filter steps")

        start_time = time.time()
        query = f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta"
        logger.debug(
            f"SimilaritySearch ({mode}, translation={translation_method}) vs host "
            f"baits: {query} vs {hostgenesbaits}"
        )
        SimilaritySearch(query, hostgenesbaits, threads, mode, translation_method)
        logger.debug(
            f"FilterTable HOST: {query}.blastx (range_junction={range_junction})"
        )
        FilterTable(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx",
            range_junction,
            "HOST",
            outdir,
        )
        logger.debug("CompareResults: dropping EEs that hit host baits harder")
        CompareResults(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred",
        )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Filter step",
                start_time=start_time,
                end_time=end_time,
                message=f"Filter step based on similarity analysis with {mode} was performed using {query} against {database}. "
                + f"Matches against same subject sequence in a {range_junction}nt range junction are filtered, mantaining the one with the greatest bitscore. "
                + "The results are compared, and the putative EEs with the greatest bitscore on host genes baits database are removed",
            )
        )
    except Exception as err:
        click.secho(f"Failed to filter EEs: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Getting basic taxonomy info")
        start_time = time.time()

        logger.debug(f"GetTaxonomy: joining hits to metadata {dbmetadata}")
        GetTaxonomy(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr",
            dbmetadata,
        )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Get Basic Taxonomy",
                start_time=start_time,
                end_time=end_time,
                message=f"Performed initial taxonomy",
            )
        )
    except Exception as err:
        click.secho(f"Failed to perform taxonomy signature: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Merging truncated elements")
        start_time = time.time()

        logger.debug(f"Merging truncated elements by {merge_level} within {limit} nt")
        GetAnnotBed(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax",
            merge_level,
        )
        MergeBed(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed",
            str(limit),
        )
        RemoveAnnotation(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge"
        )
        GetFasta(
            f"{outdir}/{prefix}.rn.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
        )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Merge truncated elements",
                start_time=start_time,
                end_time=end_time,
                message=f"Merge EEs near of {str(limit)}nt based on {merge_level} taxonomy information.",
            )
        )
    except Exception as err:
        click.secho(f"Failed to merge truncated elements: {err}", err=True, fg="red")
        sys.exit(1)
    if clean_masked:
        try:
            logger.info("Cleaning elements")
            start_time = time.time()

            logger.debug(f"MaskClean: removing EEs with > {mask_per}% soft-masking")
            MaskClean(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
                mask_per,
            )

            end_time = time.time()
            steps_infos.append(
                StepInfo.from_times(
                    step="Clean EEs",
                    start_time=start_time,
                    end_time=end_time,
                    message=f"EEs with {mask_per} percent of lower-case letters are removed.",
                )
            )
        except Exception as err:
            click.secho(
                f"Failed to remove EEs from soft-masked regions: {err}",
                err=True,
                fg="red",
            )
            sys.exit(1)

    try:
        logger.info("Creating final taxonomy")
        start_time = time.time()

        logger.debug("GetFinalTaxonomy + TagElements (Average_pident, overlap tags)")
        GetFinalTaxonomy(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax",
        )
        TagElements(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
        )
        if clean_masked:
            GetCleanedTaxonomy(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl",
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
            )
            TagElements(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.tax",
            )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Create Final Taxonomy",
                start_time=start_time,
                end_time=end_time,
                message=f"Performed final taxonomy.",
            )
        )
    except Exception as err:
        click.secho(f"Failed to obtain final taxonomy: {err}", err=True, fg="red")
        sys.exit(1)

    if overlap != "keep":
        try:
            logger.info(f"Filtering overlaping elements ({overlap})")
            start_time = time.time()
            tmp_outputs = f"{outdir}/tmp_outputs"
            os.makedirs(tmp_outputs, exist_ok=True)

            logger.debug(
                f"FilterOverlap strategy={overlap} "
                f"target_families={list(target_families)} "
                f"non_target_families={list(non_target_families)}; "
                f"removed elements -> {tmp_outputs}"
            )
            FilterOverlap(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
                overlap,
                list(target_families),
                f"{tmp_outputs}/{prefix}.EEs.removed.fa",
                f"{tmp_outputs}/{prefix}.EEs.removed.tax.tsv",
                non_target_families=list(non_target_families),
            )
            if clean_masked:
                FilterOverlap(
                    f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl",
                    f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.tax",
                    overlap,
                    list(target_families),
                    f"{tmp_outputs}/{prefix}.EEs.cleaned.removed.fa",
                    f"{tmp_outputs}/{prefix}.EEs.cleaned.removed.tax.tsv",
                    non_target_families=list(non_target_families),
                )

            end_time = time.time()
            steps_infos.append(
                StepInfo.from_times(
                    step="Filter overlaping elements",
                    start_time=start_time,
                    end_time=end_time,
                    message=f"Resolved overlaping elements with the '{overlap}' "
                    f"strategy; filtered-out elements saved to {tmp_outputs}.",
                )
            )
        except Exception as err:
            click.secho(
                f"Failed to filter overlaping elements: {err}", err=True, fg="red"
            )
            sys.exit(1)

    try:
        logger.info("Generating GFF3 annotation")
        start_time = time.time()

        logger.debug(f"WriteGFF3 (analysis={analysis}) for prefix {prefix!r}")
        WriteGFF3(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.gff3",
            prefix=prefix,
            analysis=analysis,
        )
        if clean_masked:
            WriteGFF3(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.tax",
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.gff3",
                prefix=prefix,
                analysis=analysis,
            )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Generate GFF3 annotation",
                start_time=start_time,
                end_time=end_time,
                message="Generated GFF3 annotation of endogenous elements.",
            )
        )
    except Exception as err:
        click.secho(f"Failed to generate GFF3 annotation: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Extracting flaking regions")
        start_time = time.time()

        logger.debug(f"Extracting {flank} nt flanks around each EE")
        GetLength(f"{outdir}/{prefix}.rn.fmt")
        GetBed(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
        )
        BedFlank(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed",
            f"{outdir}/{prefix}.rn.fmt.rn.fmt.lenght",
            flank,
        )
        GetFasta(
            f"{outdir}/{prefix}.rn.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.fasta",
        )

        end_time = time.time()
        steps_infos.append(
            StepInfo.from_times(
                step="Extract flanking regions",
                start_time=start_time,
                end_time=end_time,
                message=f"Extracted {flank}nt of each flanking region of EEs.",
            )
        )

    except Exception as err:
        click.secho(f"Failed to extract flanking regions: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Organizing final outputs")
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
            f"{outdir}/{prefix}.EEs.fa",
        )
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
            f"{outdir}/{prefix}.EEs.tax.tsv",
        )
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.gff3",
            f"{outdir}/{prefix}.EEs.gff3",
        )
        if clean_masked:
            os.rename(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl",
                f"{outdir}/{prefix}.EEs.cleaned.fa",
            )
            os.rename(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.tax",
                f"{outdir}/{prefix}.EEs.cleaned.tax.tsv",
            )
            os.rename(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.gff3",
                f"{outdir}/{prefix}.EEs.cleaned.gff3",
            )
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.fasta",
            f"{outdir}/{prefix}.EEs.flanks.fa",
        )
        print("")
        print("Output files:\n")
        print(
            f"{outdir}/{prefix}.EEs.fa ----------------------------- Fasta file with Endogenous Elements nucleotide sequences."
        )
        print(
            f"{outdir}/{prefix}.EEs.tax.tsv ------------------------ TSV file with Endogenous Elements taxonomy."
        )
        print(
            f"{outdir}/{prefix}.EEs.gff3 --------------------------- GFF3 annotation of Endogenous Elements."
        )
        print(
            f"{outdir}/{prefix}.EEs.flanks.fa ---------------------- Fasta file with Endogenous Elements plus {flank}nt in each flanking regions."
        )
        if clean_masked:
            print(
                f"{outdir}/{prefix}.EEs.cleaned.fa --------------------- Fasta file with Cleaned Endogenous Elements."
            )
            print(
                f"{outdir}/{prefix}.EEs.cleaned.tax.tsv ---------------- TSV file with Cleaned Endogenous Elements."
            )
            print(
                f"{outdir}/{prefix}.EEs.cleaned.gff3 ------------------- GFF3 annotation of Cleaned Endogenous Elements."
            )
        print("")
        if removetmp:
            logger.warning("Removing temporary files.\n")
            for tmp_file in glob.glob(f"{outdir}/*rn*"):
                logger.debug(f"Removing temporary file {tmp_file}")
                os.remove(tmp_file)
        else:
            if os.path.isdir(f"{outdir}/tmp_files") == False:
                os.mkdir(f"{outdir}/tmp_files")
            else:
                pass
            for tmp_file in glob.glob(f"{outdir}/*rn*"):
                new_tmp_file = re.sub(r".*/", "", tmp_file)
                os.rename(tmp_file, f"{outdir}/tmp_files/{new_tmp_file}")
            logger.info(
                f"Temporary files were moved to {outdir}/tmp_files. Check the tool documentation to access the description of each temporary file.\n"
            )
        print_info.print_finish()
    except Exception as err:
        click.secho(f"Failed to organize outputs: {err}", err=True, fg="red")
        sys.exit(1)
    end_running_time = time.time()
    run_arguments = RunArguments(
        genome_file=genome_file,
        prefix=prefix,
        outdir=outdir,
        database=database,
        dbmetadata=dbmetadata,
        baits=hostgenesbaits,
        mode=mode,
        length=length,
        flank=flank,
        limit=limit,
        range_junction=range_junction,
        mask_per=mask_per,
        clean_masked=clean_masked,
        threads=threads,
        removetmp=removetmp,
        index_databases=index_databases,
        merge_level=merge_level,
        analysis=analysis,
        overlap=overlap,
        target_families=list(target_families),
        non_target_families=list(non_target_families),
        translation_method=translation_method,
    )
    run_info = RunInfo.from_run(
        __version__,
        system_info,
        run_arguments,
        dependencies,
        start_running_time,
        end_running_time,
        steps_infos,
    )
    logger.debug(f"Writing run summary to {outdir}/eefinder.log")
    with open(f"{outdir}/eefinder.log", "w") as json_out:
        json.dump(asdict(run_info), json_out, indent=4)


def _common_download_options(func):
    """Attach the -od/-pr/--refseq/--debug options shared by every command."""
    func = click.option(
        "--stall-timeout",
        help="Seconds without any sign of life (no output from datasets, no "
        "growth of the archive) after which a download is treated as hung, "
        f"killed and retried. 0 waits forever. default = {DEFAULT_STALL_TIMEOUT}",
        type=float,
        default=DEFAULT_STALL_TIMEOUT,
    )(func)
    func = click.option(
        "--released-before",
        help="Only include data released on or before this date (YYYY-MM-DD), so "
        "a database can be rebuilt as it stood then. For bacteria/host the NCBI "
        "client applies it per assembly; for viruses its subcommand has no such "
        "flag, so EEfinder applies it per organism, from the release dates in the "
        "download's own report.",
        default=None,
    )(func)
    func = click.option(
        "--keep-download/--remove-download",
        help="Keep the downloaded zip and the extracted ncbi_dataset directory. "
        "They are deleted by default, since everything in them is already in the "
        "FASTA, the metadata CSV and the tracking table. default = remove",
        default=False,
    )(func)
    func = click.option(
        "--attempts",
        help="How many times to try the download before giving up. NCBI "
        f"transfers fail and hang intermittently. default = {DEFAULT_ATTEMPTS}",
        type=int,
        default=DEFAULT_ATTEMPTS,
    )(func)
    func = click.option(
        "--debug",
        help="Emit verbose debug logging (download command, extraction, "
        "per-record standardization details). default = off",
        is_flag=True,
    )(func)
    func = click.option(
        "--cluster/--no-cluster",
        help="Collapse 100%-identical / 100%-coverage duplicate proteins with "
        "cd-hit before writing the database. default = cluster",
        default=True,
    )(func)
    func = click.option(
        "--exclude-taxon",
        "exclude_taxa",
        help="Tax id or scientific name to leave out of the download; repeatable. "
        "The branch is never requested from NCBI: the taxon being downloaded is "
        "expanded into the taxa that cover it minus this one. 'virus' defaults "
        f"to {DEFAULT_VIRUS_EXCLUSIONS[0]} (SARS-CoV-2), which is 61% of every "
        f"viral record in GenBank; pass '--exclude-taxon {NO_EXCLUSION}' to "
        "switch that off.",
        multiple=True,
    )(func)
    func = click.option(
        "--refseq/--all-sequences",
        help="Restrict the download to RefSeq (default) or fetch all sequences.",
        default=True,
    )(func)
    func = click.option(
        "-pr",
        "--prefix",
        help="Basename for the output files, default = the dataset type "
        "(e.g. virus.fa/virus.csv).",
        default=None,
    )(func)
    func = click.option(
        "-od",
        "--outdir",
        help="Path and dir to store the downloaded database.",
        required=True,
    )(func)
    return func


def _resolve_exclusions(exclude_taxa, defaults=()):
    """Apply the per-dataset default exclusions and the 'none' escape hatch.

    Passing ``--exclude-taxon none`` clears the defaults, so a build that really
    wants every branch can still ask for one.
    """
    given = tuple(t.strip() for t in (exclude_taxa or ()) if t.strip())
    if not given:
        return tuple(defaults)
    if any(t.lower() == NO_EXCLUSION for t in given):
        kept = tuple(t for t in given if t.lower() != NO_EXCLUSION)
        return kept
    return given


def _run_get_databases(
    dataset,
    taxon,
    outdir,
    prefix,
    refseq,
    exclude_uninformative,
    standardize_proteins,
    cluster=True,
    debug=False,
    attempts=DEFAULT_ATTEMPTS,
    stall_timeout=DEFAULT_STALL_TIMEOUT,
    keep_download=False,
    released_before=None,
    exclude_taxa=(),
):
    """Check for the datasets binary and run :class:`GetDatabases`."""
    if debug:
        enable_debug()
    logger.debug(
        f"get-databases {dataset} arguments: taxon={taxon!r} outdir={outdir!r} "
        f"prefix={prefix!r} refseq={refseq} "
        f"exclude_uninformative={exclude_uninformative} "
        f"standardize_proteins={standardize_proteins} cluster={cluster} "
        f"attempts={attempts} stall_timeout={stall_timeout} "
        f"keep_download={keep_download} released_before={released_before!r} "
        f"exclude_taxa={exclude_taxa!r}"
    )
    if shutil.which(DATASETS_BINARY) is None:
        click.secho(
            f"'{DATASETS_BINARY}' was not found on PATH. Install the NCBI "
            "datasets CLI (conda package 'ncbi-datasets-cli', pinned in env.yml).",
            err=True,
            fg="red",
        )
        sys.exit(1)
    if cluster and shutil.which(CDHIT_BINARY) is None:
        click.secho(
            f"'{CDHIT_BINARY}' was not found on PATH. Install it (conda package "
            "'cd-hit', pinned in env.yml) or pass --no-cluster.",
            err=True,
            fg="red",
        )
        sys.exit(1)
    try:
        GetDatabases(
            dataset=dataset,
            taxon=taxon,
            outdir=outdir,
            prefix=prefix,
            refseq=refseq,
            exclude_uninformative=exclude_uninformative,
            standardize_proteins=standardize_proteins,
            cluster=cluster,
            attempts=attempts,
            stall_timeout=stall_timeout or None,
            keep_download=keep_download,
            released_before=released_before,
            exclude_taxa=exclude_taxa,
        )
    except Exception as err:
        click.secho(f"Failed to download databases: {err}", err=True, fg="red")
        sys.exit(1)


@cli.group(name="get-databases")
@click.version_option(__version__)
def get_databases():
    """Download RefSeq protein databases (and metadata) via NCBI datasets.

    Each group has its own subcommand: 'virus' and 'bacteria' produce a protein
    FASTA + metadata CSV (the screening -db/-mt inputs); 'host' produces the -bt
    baits FASTA (no CSV).
    """


@get_databases.command(name="virus")
@_common_download_options
@click.option(
    "-tx",
    "--taxon",
    help="NCBI taxon name or tax id to download (e.g. Flaviviridae, 10239). "
    "default = 10239 (Viruses).",
    default=10239,
    type=str,
)
@click.option(
    "--exclude-uninformative/--keep-uninformative",
    help="Drop 'hypothetical protein' and 'uncharacterized protein' records "
    "from the downloaded database. default = exclude",
    default=True,
)
@click.option(
    "--standardize-proteins/--raw-proteins",
    help="Rewrite the metadata CSV 'Protein' column to canonical names using "
    "the bundled viral protein map (also removes special characters and "
    "capitalises the first letter). default = standardize",
    default=True,
)
def get_databases_virus(
    outdir,
    prefix,
    cluster,
    refseq,
    attempts,
    stall_timeout,
    keep_download,
    released_before,
    debug,
    taxon,
    exclude_uninformative,
    standardize_proteins,
    exclude_taxa,
):
    """Download the RefSeq viral protein DB + metadata CSV (screening -db/-mt)."""
    _run_get_databases(
        dataset="virus",
        taxon=taxon or DEFAULT_TAXA["virus"],
        outdir=outdir,
        prefix=prefix,
        refseq=refseq,
        exclude_uninformative=exclude_uninformative,
        standardize_proteins=standardize_proteins,
        cluster=cluster,
        debug=debug,
        attempts=attempts,
        stall_timeout=stall_timeout,
        keep_download=keep_download,
        released_before=released_before,
        exclude_taxa=_resolve_exclusions(exclude_taxa, DEFAULT_VIRUS_EXCLUSIONS),
    )


@get_databases.command(name="bacteria")
@_common_download_options
@click.option(
    "-tx",
    "--taxon",
    help="NCBI taxon name or tax id to download (e.g. Rickettsiales, 2). "
    "default = 2 (Bacteria).",
    default=2,
    type=str,
)
@click.option(
    "--exclude-uninformative/--keep-uninformative",
    help="Drop 'hypothetical protein' and 'uncharacterized protein' records "
    "from the downloaded database. default = exclude",
    default=True,
)
@click.option(
    "--standardize-proteins/--raw-proteins",
    help="Clean the metadata CSV 'Protein' column (remove NCBI '[key=value]' "
    "tags and special characters, fix common misspellings, capitalise the first "
    "letter, drop bare CDS/ORF records). No bacterial name map is applied yet. "
    "default = standardize",
    default=True,
)
def get_databases_bacteria(
    outdir,
    prefix,
    cluster,
    refseq,
    attempts,
    stall_timeout,
    keep_download,
    released_before,
    debug,
    taxon,
    exclude_uninformative,
    standardize_proteins,
    exclude_taxa,
):
    """Download the RefSeq bacterial protein DB + metadata CSV (screening -db/-mt)."""
    _run_get_databases(
        dataset="bacteria",
        taxon=taxon or DEFAULT_TAXA["bacteria"],
        outdir=outdir,
        prefix=prefix,
        refseq=refseq,
        exclude_uninformative=exclude_uninformative,
        standardize_proteins=standardize_proteins,
        cluster=cluster,
        debug=debug,
        attempts=attempts,
        stall_timeout=stall_timeout,
        keep_download=keep_download,
        released_before=released_before,
        exclude_taxa=_resolve_exclusions(exclude_taxa),
    )


@get_databases.command(name="host")
@_common_download_options
@click.option(
    "-tx",
    "--taxon",
    help="NCBI taxon name or tax id of the host to download (e.g. 'Aedes "
    "aegypti', 7159). Required.",
    required=True,
    type=str,
)
@click.option(
    "--exclude-uninformative/--keep-uninformative",
    help="Drop 'hypothetical protein' and 'uncharacterized protein' records "
    "from the downloaded database. default = exclude",
    default=True,
)
def get_databases_host(
    outdir,
    prefix,
    cluster,
    refseq,
    attempts,
    stall_timeout,
    keep_download,
    released_before,
    debug,
    taxon,
    exclude_uninformative,
    exclude_taxa,
):
    """Download the host protein baits FASTA (screening -bt); no metadata CSV."""
    _run_get_databases(
        dataset="host",
        taxon=taxon,
        outdir=outdir,
        prefix=prefix,
        refseq=refseq,
        exclude_uninformative=exclude_uninformative,
        standardize_proteins=False,
        cluster=cluster,
        debug=debug,
        attempts=attempts,
        stall_timeout=stall_timeout,
        keep_download=keep_download,
        released_before=released_before,
        exclude_taxa=_resolve_exclusions(exclude_taxa),
    )
