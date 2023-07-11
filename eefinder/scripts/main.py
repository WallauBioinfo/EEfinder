#!/usr/bin/python3
# -*- coding: utf-8 -*-

import time
import click
import re
import os
import glob
import sys
import json
from eefinder.log import logger
from eefinder.run_message import PaperInfo
from eefinder.utils import check_outdir, step_info, running_info
from eefinder.prepare_data import InsertPrefix
from eefinder.clean_data import RemoveShortSequences, MaskClean
from eefinder.make_database import MakeDB
from eefinder.similarity_analysis import SimilaritySearch
from eefinder.filter_table import FilterTable
from eefinder.get_taxonomy import GetTaxonomy, GetFinalTaxonomy, GetCleanedTaxonomy
from eefinder.bed import GetFasta, GetAnnotBed, RemoveAnnotation, MergeBed, BedFlank, GetBed
from eefinder.compare_results import CompareResults
from eefinder.get_length import GetLength
from eefinder.tag_elements import TagElements
from eefinder import __version__


@click.group()
def cli():
    "This tool predict regions of Endogenous Elements in Eukaryote Genomes."
    pass


@cli.command()
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
    "--baits",
    help="Bait proteins, used to filter putative EEs .fasta file.",
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
    help="Taxonomy level to merge elements by genus or family, default = genus",
    default="genus",
    type=click.Choice(["family", "genus"]),
)
def main(
    genome_file,
    outdir,
    database,
    dbmetadata,
    baits,
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
):
    """
    This tool predict regions of Endogenous Elements in Eukaryote Genomes.
    """
    steps_infos = []
    start_running_time = time.time()
    print_info = PaperInfo()
    print_info.print_start(__version__)

    if prefix is None:
        try:
            logger.info(f"Creating prefix")
            prefix = genome_file
            prefix = re.sub("\..*", "", prefix)
            prefix = re.sub(".*/", "", prefix).rstrip("\n")
            prefix = re.sub(".*/", "", prefix).rstrip("\n")
        except Exception as err:
            click.secho(f"Failed to create prefix: {err}", err=True, fg="red")
            sys.exit(1)

    try:
        logger.info(f"Creating output directory")
        outdir = check_outdir(outdir)

    except Exception as err:
        click.secho(f"Failed to create output dir: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info(f"Preparing input data")
        start_time = time.time()

        InsertPrefix(genome_file, prefix, outdir)
        RemoveShortSequences(f"{outdir}/{prefix}.rn", length)

        end_time = time.time()
        steps_infos.append(
            step_info(
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

            MakeDB(mode, database, "prot", threads)
            MakeDB(mode, baits, "prot", threads)

            end_time = time.time()
            steps_infos.append(
                step_info(
                    step="Index databases",
                    start_time=start_time,
                    end_time=end_time,
                    message=f"Index {database} and {baits} databases.",
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

        SimilaritySearch(query, database, threads, mode)
        FilterTable(f"{query}.blastx", range_junction, "EE", outdir)

        end_time = time.time()
        steps_infos.append(
            step_info(
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

        GetFasta(
            f"{outdir}/{prefix}.rn.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta",
        )

        end_time = time.time()
        steps_infos.append(
            step_info(
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
        SimilaritySearch(query, baits, threads, mode)
        FilterTable(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx",
            range_junction,
            "HOST",
            outdir,
        )
        CompareResults(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred",
        )

        end_time = time.time()
        steps_infos.append(
            step_info(
                step="Filter step",
                start_time=start_time,
                end_time=end_time,
                message=f"Filter step based on similarity analysis with {mode} was performed using {query} against {database}. "
                + f"Matches against same subject sequence in a {range_junction}nt range junction are filtered, mantaining the one with the greatest bitscore. "
                + "The results are compared, and the putative EEs with the greatest bitscore on baits database are removed",
            )
        )
    except Exception as err:
        click.secho(f"Failed to filter EEs: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Getting basic taxonomy info")
        start_time = time.time()

        GetTaxonomy(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr",
            dbmetadata,
        )

        end_time = time.time()
        steps_infos.append(
            step_info(
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
            step_info(
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

            MaskClean(
                f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
                mask_per,
            )

            end_time = time.time()
            steps_infos.append(
                step_info(
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
            step_info(
                step="Create Final Taxonomy",
                start_time=start_time,
                end_time=end_time,
                message=f"Performed final taxonomy.",
            )
        )
    except Exception as err:
        click.secho(f"Failed to obtain final taxonomy: {err}", err=True, fg="red")
        sys.exit(1)

    try:
        logger.info("Extracting flaking regions")
        start_time = time.time()

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
            step_info(
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
            f"{outdir}/{prefix}.EEs.flanks.fa ---------------------- Fasta file with Endogenous Elements plus {flank}nt in each flanking regions."
        )
        if clean_masked:
            print(
                f"{outdir}/{prefix}.EEs.cleaned.fa --------------------- Fasta file with Cleaned Endogenous Elements."
            )
            print(
                f"{outdir}/{prefix}.EEs.cleaned.tax.tsv ---------------- TSV file with Cleaned Endogenous Elements."
            )
        print("")
        if removetmp:
            logger.warning("Removing temporary files.\n")
            for tmp_file in glob.glob(f"{outdir}/*rn*"):
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
    arguments = [
        genome_file,
        outdir,
        database,
        dbmetadata,
        baits,
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
    ]
    running_infos = running_info(
        arguments, start_running_time, end_running_time, steps_infos
    )
    with open(f"{outdir}/eefinder.log", "w") as json_out:
        json.dump(running_infos, json_out, indent=4)
