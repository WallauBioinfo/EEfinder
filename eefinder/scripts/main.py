#!/usr/bin/python3
# -*- coding: utf-8 -*-

import time
import click
import re
import os
import glob
import sys
from eefinder.run_message import PaperInfo
from eefinder.prepare_data import SetPrefix
from eefinder.clean_data import RemoveShort, MaskClean
from eefinder.make_database import MakeDB
from eefinder.similarity_analysis import SimilaritySearch
from eefinder.filter_table import FilterTable
from eefinder.get_fasta import GetFasta
from eefinder.get_taxonomy import GetTaxonomy, GetFinalTaxonomy, GetCleanedTaxonomy
from eefinder.format_bed import GetAnnotBed, RemoveAnnotation
from eefinder.bed_merge import MergeBed
from eefinder.bed_flank import BedFlank
from eefinder.compare_results import CompareResults
from eefinder.get_bed import GetBed
from eefinder.get_length import GetLength
from eefinder.tag_elements import TagElements
from eefinder import __version__


@click.group()
def cli():
    "This tool predict regions of Endogenous Elements in Eukaryote Genomes."
    pass


@cli.command()
@click.version_option(__version__)
@click.option("-in", "--input_file", help="Fasta file (nucleotides).", required=True)
@click.option("-od", "--outdir", help="Path and dir to store output", required=True)
@click.option(
    "-db",
    "--database",
    help="Proteins from viruses or bacterias database fasta file.",
    required=True,
)
@click.option(
    "-mt",
    "--dbmetadata",
    help="Proteins from viruses or bacterias csv file.",
    required=True,
)
@click.option(
    "-dbh",
    "--dbhost",
    help="Host and heatshock proteins database file.",
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
    help="Minimum length of contigs used for BLAST or DIAMOND, default = 10,000.",
    type=int,
    default=10000,
)
@click.option(
    "-fl",
    "--flank",
    help="Length of flanking regions of Endogenous Elements to be extracted, default = 10,000.",
    type=int,
    default=10000,
)
@click.option(
    "-lm",
    "--limit",
    help="Limit of bases used to merge regions on bedtools merge, default = 1.",
    type=str,
    default="1",
)
@click.option(
    "-rj,",
    "--range_junction",
    help="Sets the range for junction of redudant hits, should follow a logic with 'limit' option, default=100",
    type=int,
    default=100,
)
@click.option(
    "-mp",
    "--mask_per",
    help="Limit of lowercase letters in percentage to consider a putative Endogenous Elements as a repetitive region, default = 50.",
    type=str,
    default=50,
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
    help="Remove temporary files generated through analysis? default = True.",
    default=False,
    is_flag=True,
)
@click.option(
    "-mk",
    "--makeblastdb",
    help="Make blast database?, default = True.",
    default=False,
    is_flag=True,
)
@click.option(
    "-pr",
    "--prefix",
    help="Write the prefix name for output files. We strongly recommend the use of genome name assembly, this prefix will be used to create the EEname (The Endogenous Element name will be formated as PREFIX|CONTIG/SCAFFOLD:START-END) default = input file name.",
)
@click.option(
    "-ml",
    "--merge_level",
    help="Phylogenetic level to merge elements by genus or family, default = family",
    default="genus",
    type=click.Choice(["family", "genus"]),
)
def main(
    input_file,
    outdir,
    database,
    dbmetadata,
    dbhost,
    mode,
    length,
    flank,
    limit,
    range_junction,
    mask_per,
    threads,
    removetmp,
    makeblastdb,
    prefix,
    merge_level,
):  
    if prefix is not None:
        try:
            prefix = input_file
            prefix = re.sub("\..*", "", prefix)
            prefix = re.sub(".*/", "", prefix).rstrip("\n")
            prefix = re.sub(".*/", "", prefix).rstrip("\n")
        except Exception as err:
            click.secho(f"Failed to create prefix: {err}", err=True, fg="red")
            sys.exit(1)
    try:
        if outdir != None:
            if "/" in outdir:
                outdir = re.sub("/$", "", outdir)
        else:
            if "/" not in input_file:
                outdir = "."
            else:
                outdir = re.sub("/.*", "", input_file).rstrip("\n")
        if os.path.isdir(outdir) == False:
            os.mkdir(outdir)
    except Exception as err:
        click.secho(f"Failed to create output dir: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        log_file = open(f"{outdir}/EEfinder.log.txt", "w+")
        start_running_time = time.time()
        print_info = PaperInfo()
        print_info.print_message()
        # Prepare Data Step
        print("|" + "-" * 45 + "PREPARING DATA" + "-" * 45 + "|\n")
        print("|" + "-" * 45 + "PREPARING DATA" + "-" * 45 + "|\n", file=log_file)
        start_time_prep = time.time()

        # Set Prefix
        set_prefix = SetPrefix(input_file, prefix, outdir)
        set_prefix.run_insert_prefix()
        # Remove short sequences
        remove_seqs = RemoveShort(f"{outdir}/{prefix}.rn", length, log_file)
        remove_seqs.run_remove_short()
        final_time_prep = time.time()
        print(
            f"PREPARING DATA TIME = {(final_time_prep - start_time_prep)/60:.2f} MINUTES"
        )
        print(
            f"PREPARING DATA TIME = {(final_time_prep - start_time_prep)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to prepare input data: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        print("\n|" + "-" * 42 + "FORMATTING DATABASES" + "-" * 42 + "|\n")
        print(
            "\n|" + "-" * 42 + "FORMATTING DATABASES" + "-" * 42 + "|\n", file=log_file
        )
        start_time_formatdb = time.time()
        makeblastdb_ee = MakeDB(mode, database, "prot", threads, makeblastdb, log_file)
        makeblastdb_ee.run_make_db()
        makeblastdb_filter = MakeDB(
            mode, dbhost, "prot", threads, makeblastdb, log_file
        )
        makeblastdb_filter.run_make_db()
        final_time_formatdb = time.time()
        print(
            f"FORMATTING DATABASES TIME = {(final_time_formatdb - start_time_formatdb)/60:.2f} MINUTES"
        )
        print(
            f"FORMATTING DATABASES TIME = {(final_time_formatdb - start_time_formatdb)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to format databases: {err}", err=True, fg="red")
        sys.exit(1)    
    try:
        print("\n|" + "-" * 40 + "RUNNING SIMILARITY SEARCH" + "-" * 39 + "|\n")
        print(
            "\n|" + "-" * 40 + "RUNNING SIMILARITY SEARCH" + "-" * 39 + "|\n",
            file=log_file,
        )
        start_time_sim = time.time()
        ee_similarity_step = SimilaritySearch(
            f"{outdir}/{prefix}.rn.fmt", database, threads, mode, log_file
        )
        ee_similarity_step.run_similarity_search()
        # Filtering results
        ee_filter_table = FilterTable(
            f"{outdir}/{prefix}.rn.fmt.blastx", range_junction, "EE", outdir, log_file
        )
        ee_filter_table.run_filter()
        final_time_sim = time.time()
        print(
            f"RUNING SEARCH TIME = {(final_time_sim - start_time_sim)/60:.2f} MINUTES"
        )
        print(
            f"RUNING SEARCH TIME = {(final_time_sim - start_time_sim)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to run similarity searches: {err}", err=True, fg="red")
        sys.exit(1)  
    try:
        print("\n|" + "-" * 40 + "EXTRACTING PUTATIVE EVES" + "-" * 40 + "|\n")
        print(
            "\n|" + "-" * 40 + "EXTRACTING PUTATIVE EVES" + "-" * 40 + "|\n",
            file=log_file,
        )
        start_time_extract = time.time()
        getter_fasta = GetFasta(
            f"{outdir}/{prefix}.rn.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta",
            log_file,
        )
        getter_fasta.run_get_fasta()
        final_time_extract = time.time()
        print(
            f"EXTRACTING EVES TIME = {(final_time_extract - start_time_extract)/60:.2f} MINUTES"
        )
        print(
            f"EXTRACTING EVES TIME = {(final_time_extract - start_time_extract)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to run extract EE sequences: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        # Second Similarity Screen
        print("\n|" + "-" * 42 + "RUNNING FILTER STEPS" + "-" * 42 + "|\n")
        print(
            "\n|" + "-" * 42 + "RUNNING FILTER STEPS" + "-" * 42 + "|\n", file=log_file
        )
        start_time_filter = time.time()
        host_similarity_step = SimilaritySearch(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta",
            dbhost,
            threads,
            mode,
            log_file,
        )
        host_similarity_step.run_similarity_search()
        # Filtering results
        host_filter_table = FilterTable(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx",
            range_junction,
            "HOST",
            outdir,
            log_file,
        )
        host_filter_table.run_filter()
        # Comparing results
        comparer = CompareResults(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred",
            log_file,
        )
        comparer.run_comparation()
        final_time_filter = time.time()
        print(
            f"FILTER STEP TIME = {(final_time_filter - start_time_filter)/60:.2f} MINUTES"
        )
        print(
            f"FILTER STEP TIME = {(final_time_filter - start_time_filter)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to filter EEs: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        print("\n|" + "-" * 39 + "GETTING BASIC TAXONOMY INFO" + "-" * 38 + "|\n")
        print(
            "\n|" + "-" * 39 + "GETTING BASIC TAXONOMY INFO" + "-" * 38 + "|\n",
            file=log_file,
        )
        start_time_tax = time.time()
        # Get info taxonomy for pEVEs
        get_info = GetTaxonomy(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr",
            dbmetadata,
            log_file,
        )
        get_info.run_get_taxonomy_info()
        final_time_tax = time.time()
        print(
            f"GETTING BASIC TAXONOMY INFO TIME = {(final_time_tax - start_time_tax)/60:.2f} MINUTES"
        )
        print(
            f"GETTING BASIC TAXONOMY INFO TIME = {(final_time_tax - start_time_tax)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to perform taxonomy signature: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        print("\n|" + "-" * 40 + "MERGIN TRUNCATED ELEMENTS" + "-" * 39 + "|\n")
        print(
            "\n|" + "-" * 40 + "MERGIN TRUNCATED ELEMENTS" + "-" * 39 + "|\n",
            file=log_file,
        )
        start_time_merge = time.time()
        get_annot_bed = GetAnnotBed(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax",
            merge_level,
            log_file,
        )
        get_annot_bed.run_get_annotated_bed()
        # Merge bed file
        merge_bed = MergeBed(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed",
            str(limit),
            log_file,
        )
        merge_bed.run_merge_bedfile()
        # Format bed names to get_fasta step and to create a taxonomy table with merged elements
        remove_annot_bed = RemoveAnnotation(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge",
            log_file,
        )
        remove_annot_bed.run_reformat_bed()
        # Getting merged elements
        getter_merged_fasta = GetFasta(
            f"{outdir}/{prefix}.rn.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
            log_file,
        )
        getter_merged_fasta.run_get_fasta()
        final_time_merge = time.time()
        print(f"MERGING TIME = {(final_time_merge - start_time_merge)/60:.2f} MINUTES")
        print(
            f"MERGING TIME = {(final_time_merge - start_time_merge)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to merge truncated elements: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        print("\n|" + "-" * 43 + "CLEANNING ELEMENTS" + "-" * 43 + "|\n")
        print("\n|" + "-" * 43 + "CLEANNING ELEMENTS" + "-" * 43 + "|\n", file=log_file)
        start_time_clean = time.time()
        clean_masked = MaskClean(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
            mask_per,
            log_file,
        )
        clean_masked.run_mask_clean()
        final_time_clean = time.time()
        print(
            f"CLEAN STEP TIME = {(final_time_clean - start_time_clean)/60:.2f} MINUTES"
        )
        print(
            f"CLEAN STEP TIME = {(final_time_clean - start_time_clean)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to remove EEs from soft-masked regions: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        print("\n|" + "-" * 37 + "CREATING FINAL TAXONOMY FILES" + "-" * 38 + "|\n")
        print(
            "\n|" + "-" * 37 + "CREATING FINAL TAXONOMY FILES" + "-" * 38 + "|\n",
            file=log_file,
        )
        start_time_final_tax = time.time()
        # Create final taxonomy files to '.tax.bed.merge.fmt.fa' and '.tax.bed.merge.fmt.fa.mask_clean.fa'
        get_final_taxonomy = GetFinalTaxonomy(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax",
            log_file,
        )
        get_final_taxonomy.run_get_final_taxonomy()
        get_cleaned_taxonomy = GetCleanedTaxonomy(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
            log_file,
        )
        get_cleaned_taxonomy.run_get_cleaned_taxonomy()
        tag_taxonomy = TagElements(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
            log_file,
        )
        tag_taxonomy.run_tag_elemets()
        tag_cleaned_taxonomy = TagElements(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.tax",
            log_file,
        )
        tag_cleaned_taxonomy.run_tag_elemets()

        final_time_final_tax = time.time()
        print(
            f"FINAL TAXONOMY TIME = {(final_time_final_tax - start_time_final_tax)/60:.2f} MINUTES"
        )
        print(
            f"FINAL TAXONOMY TIME = {(final_time_final_tax - start_time_final_tax)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to obtain final taxonomy: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        # Flanking Regions
        print("\n|" + "-" * 39 + "EXTRACTING FLANKING REGIONS" + "-" * 38 + "|\n")
        print(
            "\n|" + "-" * 39 + "EXTRACTING FLANKING REGIONS" + "-" * 38 + "|\n",
            file=log_file,
        )
        start_time_flank = time.time()
        getter_length = GetLength(f"{outdir}/{prefix}.rn.fmt", log_file)
        getter_length.run_get_length()
        getter_bed = GetBed(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
            log_file,
        )
        getter_bed.run_get_bed()
        getter_bed_flank = BedFlank(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed",
            f"{outdir}/{prefix}.rn.fmt.rn.fmt.lenght",
            flank,
            log_file,
        )
        getter_bed_flank.run_bed_flank()
        get_fasta_flank = GetFasta(
            f"{outdir}/{prefix}.rn.fmt",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank",
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.fasta",
            log_file,
        )
        get_fasta_flank.run_get_fasta()
        final_time_flank = time.time()
        print(
            f"EXTRACTING FLANKS TIME = {(final_time_flank - start_time_flank)/60:.2f} MINUTES"
        )
        print(
            f"EXTRACTING FLANKS TIME = {(final_time_flank - start_time_flank)/60:.2f} MINUTES",
            file=log_file,
        )
    except Exception as err:
        click.secho(f"Failed to extract flanking regions: {err}", err=True, fg="red")
        sys.exit(1)
    try:
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
            f"{outdir}/{prefix}.EEs.fa",
        )
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl",
            f"{outdir}/{prefix}.EEs.cleaned.fa",
        )
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
            f"{outdir}/{prefix}.EEs.tax.tsv",
        )
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.tax",
            f"{outdir}/{prefix}.EEs.cleaned.tax.tsv",
        )
        os.rename(
            f"{outdir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.fasta",
            f"{outdir}/{prefix}.EEs.flanks.fa",
        )

        print("\n|" + "-" * 45 + "SUMMARY RESULTS" + "-" * 45 + "|\n")
        print("\n|" + "-" * 45 + "SUMMARY RESULTS" + "-" * 45 + "|\n", file=log_file)
        print(
            f"{outdir}/{prefix}.EEs.fa ----------------------------- Fasta file with Endogenous Elements nucleotide sequences."
        )
        print(
            f"{outdir}/{prefix}.EEs.tax.tsv ------------------------ TSV file with Endogenous Elements taxonomy."
        )
        print(
            f"{outdir}/{prefix}.EEs.flanks.fa ---------------------- Fasta file with Endogenous Elements plus {flank}nt in each flanking regions."
        )
        print(
            f"{outdir}/{prefix}.EEs.cleaned.fa --------------------- Fasta file with Cleaned Endogenous Elements."
        )
        print(
            f"{outdir}/{prefix}.EEs.cleaned.tax.tsv ---------------- TSV file with Cleaned Endogenous Elements."
        )

        print(
            f"{outdir}/{prefix}.EEs.fa ----------------------------- Fasta file with Endogenous Elements nucleotide sequences.",
            file=log_file,
        )
        print(
            f"{outdir}/{prefix}.EEs.tax.tsv ------------------------ TSV file with Endogenous Elements taxonomy.",
            file=log_file,
        )
        print(
            f"{outdir}/{prefix}.EEs.flanks.fa ---------------------- Fasta file with Endogenous Elements plus {flank}nt in each flanking regions.",
            file=log_file,
        )
        print(
            f"{outdir}/{prefix}.EEs.cleaned.fa --------------------- Fasta file with Cleaned Endogenous Elements.",
            file=log_file,
        )
        print(
            f"{outdir}/{prefix}.EEs.cleaned.tax.tsv ---------------- TSV file with Cleaned Endogenous Elements.",
            file=log_file,
        )

        if removetmp == True:
            for tmp_file in glob.glob(f"{outdir}/*rn*"):
                os.remove(tmp_file)
            print("\nTemporary files were removed.")
            print("\nTemporary files were removed.", file=log_file)
        else:
            if os.path.isdir(f"{outdir}/tmp_files") == False:
                os.mkdir(f"{outdir}/tmp_files")
            else:
                pass
            for tmp_file in glob.glob(f"{outdir}/*rn*"):
                new_tmp_file = re.sub(r".*/", "", tmp_file)
                os.rename(tmp_file, f"{outdir}/tmp_files/{new_tmp_file}")
            print(
                f"\nTemporary files were moved to {outdir}/tmp_files. Check the github documentation to access the description of each temporary file."
            )
            print(
                f"\nTemporary files were moved to {outdir}/tmp_files. Check the github documentation to access the description of each temporary file.",
                file=log_file,
            )
        print_info.print_finish()
        end_running_time = time.time()
    except Exception as err:
        click.secho(f"Failed to organize outputs: {err}", err=True, fg="red")
        sys.exit(1)
    total_running_time = end_running_time - start_running_time
    click.secho(f"TOTAL TIME = {total_running_time/60:.2f} MINUTES", fg="green")
    print(f"TOTAL TIME = {total_running_time/60:.2f} MINUTES", file=log_file)
