#!/usr/bin/python3
# -*- coding: utf-8 -*-

##################################>IMPORT-MODULES<##################################

import time
import argparse
import re
import os
import glob
from scripts.run_message import PaperInfo
from scripts.prepare_data import ConcatFiles, SetPrefix
from scripts.clean_data import RemoveShort, MaskClean
from scripts.make_database import MakeDB
from scripts.similarity_analysis import FlankSearch, SimilaritySearch
from scripts.filter_table import FilterTable
from scripts.get_fasta import GetFasta
from scripts.get_taxonomy import GetTaxonomy, GetFinalTaxonomy, GetCleanedTaxonomy
from scripts.format_bed import GetAnnotBed, RemoveAnnotation
from scripts.bed_merge import MergeBed
from scripts.bed_flank import BedFlank
from scripts.compare_results import CompareResults
from scripts.flank_divider import FlankDivider
from scripts.get_bed import GetBed
from scripts.get_length import GetLength

###############################>GENERAL-INFORMATIONS<###############################

__authors__ = {"Names": ["Filipe Dezordi", "Yago Dias"],
               "emails": ["zimmer.filipe@gmail.com", "yag.dias@gmail"],
               "gits": ["https://github.com/dezordi", "https://github.com/yagdias"]}
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Filipe Dezordi"
__status__ = "Prototype"

#####################################>ARGUMENTS<#####################################
parser = argparse.ArgumentParser(description='This tool predict regions of Endogenous Elements in Eukaryote Genomes.',
                                 usage='\n>Running with default parameters:\n   python EEfinder.py -in <genome.fasta> -db <viral/bac_prot.fasta> -mt <viral/bac_prot.csv> -db1 <host_proteins.fasta> -db2 <host_tranposons_proteins.fasta>\n\nFor more examples, see the documentation on Github.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-in", "--input",
                    help="Fasta file (nucleotides).", required=True)
parser.add_argument("-od", "--outdir",
                    help="Path and dir to store output", required=True)
parser.add_argument("-db", "--database",
                    help="Proteins from viruses or bacterias database fasta file.", required=True)
parser.add_argument("-mt", "--dbmetadata",
                    help="Proteins from viruses or bacterias csv file.", required=True)
parser.add_argument("-db2", '--databaseHOST',
                    help="Host and heatshock proteins database file.", required=True)
parser.add_argument("-db3", "--databaseTE",
                    help="Transposon protein database file.", required=True)
parser.add_argument("-md", "--mode",
                    help="Choose between BLAST or the DIAMOND strategies (fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensisitve) to run analysis, default = blastx.", default='blastx', choices=['blastx', 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'])
parser.add_argument("-ln", "--length",
                    help="Minimum length of contigs used for BLAST or DIAMOND, default = 10,000.", type=int, default=10000)
parser.add_argument("-fl", "--flank",
                    help="Length of flanking regions of Endogenous Elements to be extracted, default = 10,000.", type=int, default=10000)
parser.add_argument("-lm", "--limit",
                    help="Limit of bases used to merge regions on bedtools merge, default = 1.", type=str, default=str(1))
parser.add_argument("-mp", "--mask_per",
                    help="Limit of lowercase letters in percentage to consider a putative Endogenous Elements as a repetitive region, default = 50.", type=str, default=50)
parser.add_argument("-p", "--threads",
                    help="Threads for multi-thread analysis, default = 1.", type=int, default=1)
parser.add_argument("-rm", "--removetmp",
                    help="Remove temporary files generated through analysis? default = True.", default=True, choices=['True', 'False'])
parser.add_argument("-mk", "--makeblastdb",
                    help="Make blast database?, default = True.", default=True, choices=['True', 'False'])
parser.add_argument("-pr", "--prefix",
                    help="Write the prefix name for output files. We strongly recommend the use of genome name assembly, this prefix will be used to create the EEname (The Endogenous Element name will be formated as PREFIX|CONTIG/SCAFFOLD:START-END) default = input file name.")
parser.add_argument("-ml", "--merge_level", help="Phylogenetic level to merge elements by genus or family, default = family",
                    default="genus", choices=['family', 'genus'])

# Store arguments on variables
args = parser.parse_args()
input_file = args.input
out_dir = args.outdir
database_file = args.database
database_table = args.dbmetadata
database_HOST = args.databaseHOST
database_TE = args.databaseTE
mode = args.mode
length_cutoff = args.length
flank_region = args.flank
limit_merge = args.limit
mask_per = args.mask_per
threads = args.threads
temp_remove = args.removetmp
make_db = args.makeblastdb
prefix = args.prefix
merge_level = args.merge_level

# Set Prefix and StorePath
if prefix == None:  # format prefix name, if the -pr argument was not used
    prefix = input_file
    prefix = re.sub("\..*", "", prefix)
    prefix = re.sub('.*/', '', prefix).rstrip('\n')
    prefix = re.sub('.*/', '', prefix).rstrip('\n')
# Create inputdir path
in_dir = re.sub("/[^/]*$", "/", database_file).rstrip('\n')
# checar se n tiver um diretorio, criar um diretorio
if out_dir != None:
    if "/" in out_dir:
        out_dir = re.sub("/$", "", out_dir)
else:
    if '/' not in input_file:
        out_dir = '.'
    else:
        out_dir = re.sub("/.*", "", input_file).rstrip('\n')
if os.path.isdir(out_dir) == False:
    os.mkdir(out_dir)
# create log file
log_file = open(f"{out_dir}/EEfinder.log.txt", "w+")

if __name__ == '__main__':
    start_running_time = time.time()
    print_info = PaperInfo()
    print_info.print_message()
    # Prepare Data Step
    print("|"+"-"*45+"PREPARING DATA"+"-"*45+"|\n")
    print("|"+"-"*45+"PREPARING DATA"+"-"*45+"|\n", file=log_file)
    start_time_prep = time.time()
    concat_files = ConcatFiles(database_HOST,
                               database_TE,
                               in_dir,
                               'DatabaseFilter.fa',
                               log_file)
    concat_files.run_concate_file()
    # Set Prefix
    set_prefix = SetPrefix(input_file,
                           prefix,
                           out_dir)
    set_prefix.run_insert_prefix()
    # Remove short sequences
    remove_seqs = RemoveShort(f"{out_dir}/{prefix}.rn",
                              length_cutoff,
                              log_file)
    remove_seqs.run_remove_short()
    final_time_prep = time.time()
    print(
        f"PREPARING DATA TIME = {(final_time_prep - start_time_prep)/60:.2f} MINUTES")
    print(
        f"PREPARING DATA TIME = {(final_time_prep - start_time_prep)/60:.2f} MINUTES", file=log_file)
    # Make databases for EVEs screening and filters
    print("\n|"+"-"*42+"FORMATTING DATABASES"+"-"*42+"|\n")
    print("\n|"+"-"*42+"FORMATTING DATABASES"+"-"*42+"|\n", file=log_file)
    start_time_formatdb = time.time()
    make_db_ee = MakeDB(mode,
                        database_file,
                        'prot',
                        threads,
                        make_db,
                        log_file)
    make_db_ee.run_make_db()
    make_db_filter = MakeDB(mode,
                            f"{in_dir}/DatabaseFilter.fa",
                            'prot',
                            threads,
                            make_db,
                            log_file)
    make_db_filter.run_make_db()
    final_time_formatdb = time.time()
    print(
        f"FORMATTING DATABASES TIME = {(final_time_formatdb - start_time_formatdb)/60:.2f} MINUTES")
    print(
        f"FORMATTING DATABASES TIME = {(final_time_formatdb - start_time_formatdb)/60:.2f} MINUTES", file=log_file)
    # Run Similarity Screenings
    print("\n|"+"-"*40+"RUNNING SIMILARITY SEARCH"+"-"*39+"|\n")
    print("\n|"+"-"*40+"RUNNING SIMILARITY SEARCH"+"-"*39+"|\n", file=log_file)
    start_time_sim = time.time()
    ee_similarity_step = SimilaritySearch(f"{out_dir}/{prefix}.rn.fmt",
                                          database_file,
                                          threads,
                                          mode,
                                          log_file)
    ee_similarity_step.run_similarity_search()
    # Filtering results
    ee_filter_table = FilterTable(f"{out_dir}/{prefix}.rn.fmt.blastx",
                                  "EE",
                                  log_file)
    ee_filter_table.run_filter()
    final_time_sim = time.time()
    print(
        f"RUNING SEARCH TIME = {(final_time_sim - start_time_sim)/60:.2f} MINUTES")
    print(
        f"RUNING SEARCH TIME = {(final_time_sim - start_time_sim)/60:.2f} MINUTES", file=log_file)
    # Getting fastas
    print("\n|"+"-"*40+"EXTRACTING PUTATIVE EVES"+"-"*40+"|\n")
    print("\n|"+"-"*40+"EXTRACTING PUTATIVE EVES"+"-"*40+"|\n", file=log_file)
    start_time_extract = time.time()
    getter_fasta = GetFasta(f"{out_dir}/{prefix}.rn.fmt",
                            f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed",
                            f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta",
                            log_file)
    getter_fasta.run_get_fasta()
    final_time_extract = time.time()
    print(
        f"EXTRACTING EVES TIME = {(final_time_extract - start_time_extract)/60:.2f} MINUTES")
    print(
        f"EXTRACTING EVES TIME = {(final_time_extract - start_time_extract)/60:.2f} MINUTES", file=log_file)
    # Second Similarity Screen
    print("\n|"+"-"*42+"RUNNING FILTER STEPS"+"-"*42+"|\n")
    print("\n|"+"-"*42+"RUNNING FILTER STEPS"+"-"*42+"|\n", file=log_file)
    start_time_filter = time.time()
    host_similarity_step = SimilaritySearch(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta",
                                            f"{in_dir}DatabaseFilter.fa",
                                            threads,
                                            mode,
                                            log_file)
    host_similarity_step.run_similarity_search()
    # Filtering results
    host_filter_table = FilterTable(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx",
                                    "HOST",
                                    log_file)
    host_filter_table.run_filter()
    # Comparing results
    comparer = CompareResults(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred",
                              f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred",
                              log_file)
    comparer.run_comparation()
    final_time_filter = time.time()
    print(
        f"FILTER STEP TIME = {(final_time_filter - start_time_filter)/60:.2f} MINUTES")
    print(
        f"FILTER STEP TIME = {(final_time_filter - start_time_filter)/60:.2f} MINUTES", file=log_file)
    print("\n|"+"-"*39+"GETTING BASIC TAXONOMY INFO"+"-"*38+"|\n")
    print("\n|"+"-"*39+"GETTING BASIC TAXONOMY INFO"+"-"*38+"|\n", file=log_file)
    start_time_tax = time.time()
    # Get info taxonomy for pEVEs
    get_info = GetTaxonomy(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr",
                           database_table,
                           log_file)
    get_info.run_get_taxonomy_info()
    final_time_tax = time.time()
    print(
        f"GETTING BASIC TAXONOMY INFO TIME = {(final_time_tax - start_time_tax)/60:.2f} MINUTES")
    print(
        f"GETTING BASIC TAXONOMY INFO TIME = {(final_time_tax - start_time_tax)/60:.2f} MINUTES", file=log_file)
    # Get annotated bed file
    print("\n|"+"-"*40+"MERGIN TRUNCATED ELEMENTS"+"-"*39+"|\n")
    print("\n|"+"-"*40+"MERGIN TRUNCATED ELEMENTS"+"-"*39+"|\n", file=log_file)
    start_time_merge = time.time()
    get_annot_bed = GetAnnotBed(
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax", merge_level, log_file)
    get_annot_bed.run_get_annotated_bed()
    # Merge bed file
    merge_bed = MergeBed(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed",
                         str(limit_merge),
                         log_file)
    merge_bed.run_merge_bedfile()
    # Format bed names to get_fasta step and to create a taxonomy table with merged elements
    remove_annot_bed = RemoveAnnotation(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge",
                                        log_file)
    remove_annot_bed.run_reformat_bed()
    # Getting merged elements
    getter_merged_fasta = GetFasta(f"{out_dir}/{prefix}.rn.fmt",
                                   f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt",
                                   f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
                                   log_file)
    getter_merged_fasta.run_get_fasta()
    final_time_merge = time.time()
    print(
        f"MERGING TIME = {(final_time_merge - start_time_merge)/60:.2f} MINUTES")
    print(
        f"MERGING TIME = {(final_time_merge - start_time_merge)/60:.2f} MINUTES", file=log_file)
    print("\n|"+"-"*43+"CLEANNING ELEMENTS"+"-"*43+"|\n")
    print("\n|"+"-"*43+"CLEANNING ELEMENTS"+"-"*43+"|\n", file=log_file)
    start_time_clean = time.time()
    clean_masked = MaskClean(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
                             mask_per,
                             log_file)
    clean_masked.run_mask_clean()
    final_time_clean = time.time()
    print(
        f"CLEAN STEP TIME = {(final_time_clean - start_time_clean)/60:.2f} MINUTES")
    print(
        f"CLEAN STEP TIME = {(final_time_clean - start_time_clean)/60:.2f} MINUTES", file=log_file)
    print("\n|"+"-"*37+"CREATING FINAL TAXONOMY FILES"+"-"*38+"|\n")
    print("\n|"+"-"*37+"CREATING FINAL TAXONOMY FILES"+"-"*38+"|\n", file=log_file)
    start_time_final_tax = time.time()
    # Create final taxonomy files to '.tax.bed.merge.fmt.fa' and '.tax.bed.merge.fmt.fa.mask_clean.fa'
    get_final_taxonomy = GetFinalTaxonomy(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt",
                                          f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax",
                                          log_file)
    get_final_taxonomy.run_get_final_taxonomy()
    get_cleaned_taxonomy = GetCleanedTaxonomy(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl",
                                              f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
                                              log_file)
    get_cleaned_taxonomy.run_get_cleaned_taxonomy()

    final_time_final_tax = time.time()
    print(
        f"FINAL TAXONOMY TIME = {(final_time_final_tax - start_time_final_tax)/60:.2f} MINUTES")
    print(
        f"FINAL TAXONOMY TIME = {(final_time_final_tax - start_time_final_tax)/60:.2f} MINUTES", file=log_file)
    # Flanking Regions
    print("\n|"+"-"*39+"EXTRACTING FLANKING REGIONS"+"-"*38+"|\n")
    print("\n|"+"-"*39+"EXTRACTING FLANKING REGIONS"+"-"*38+"|\n", file=log_file)
    start_time_flank = time.time()
    getter_length = GetLength(f"{out_dir}/{prefix}.rn.fmt",
                              log_file)
    getter_length.run_get_length()
    getter_bed = GetBed(f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa',
                        log_file)
    getter_bed.run_get_bed()
    getter_bed_flank = BedFlank(f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed',
                                f'{out_dir}/{prefix}.rn.fmt.rn.fmt.lenght',
                                flank_region,
                                log_file)
    getter_bed_flank.run_bed_flank()
    divider = FlankDivider(
        f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank')
    divider.run_flank_divider()
    get_fasta_flank = GetFasta(
        f"{out_dir}/{prefix}.rn.fmt",
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.fmt",
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.fmt.fasta",
        log_file)
    get_fasta_flankr = GetFasta(
        f"{out_dir}/{prefix}.rn.fmt",
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right",
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right.fasta",
        log_file)
    get_fasta_flankl = GetFasta(
        f"{out_dir}/{prefix}.rn.fmt",
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left",
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left.fasta",
        log_file)
    get_fasta_flank.run_get_fasta()
    get_fasta_flankl.run_get_fasta()
    get_fasta_flankr.run_get_fasta()
    final_time_flank = time.time()
    print(
        f"EXTRACTING FLANKS TIME = {(final_time_flank - start_time_flank)/60:.2f} MINUTES")
    print(
        f"EXTRACTING FLANKS TIME = {(final_time_flank - start_time_flank)/60:.2f} MINUTES", file=log_file)
    print("\n|"+"-"*38+"IDENTIFYING FLANKING ELEMENTS"+"-"*37+"|\n")
    print("\n|"+"-"*38+"IDENTIFYING FLANKING ELEMENTS"+"-"*37+"|\n", file=log_file)
    start_time_flanks_sim = time.time()
    # Flank transposon search
    # Creating dbs
    make_db_tel = MakeDB('tblastn',
                         f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left.fasta',
                         'nucl',
                         threads,
                         make_db,
                         log_file)
    make_db_tel.run_make_db()

    make_db_ter = MakeDB('tblastn',
                         f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right.fasta',
                         'nucl',
                         threads,
                         make_db,
                         log_file)
    make_db_ter.run_make_db()
    # Left side
    tel_similarity = FlankSearch(database_TE,
                                 f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left.fasta',
                                 threads,
                                 f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left.fasta.tblastn',
                                 log_file)
    tel_similarity.run_flank_search()

    tel_filter = FilterTable(f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left.fasta.tblastn',
                             'HOST',
                             log_file)
    tel_filter.run_filter()
    # Right side
    ter_similarity = FlankSearch(database_TE,
                                 f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right.fasta',
                                 threads,
                                 f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right.fasta.tblastn',
                                 log_file)
    ter_similarity.run_flank_search()
    ter_filter = FilterTable(f'{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right.fasta.tblastn',
                             'HOST',
                             log_file)
    ter_filter.run_filter()
    final_time_flanks_sim = time.time()
    print(
        f"IDENTIFYING FLANKING ELEMENTS TIME = {(final_time_flanks_sim - start_time_flanks_sim)/60:.2f} MINUTES")
    print(
        f"IDENTIFYING FLANKING ELEMENTS TIME = {(final_time_flanks_sim - start_time_flanks_sim)/60:.2f} MINUTES", file=log_file)
    # Organize Outputs
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa",
              f"{out_dir}/{prefix}.EEs.fa")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl",
              f"{out_dir}/{prefix}.EEs.cleaned.fa")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.tax",
              f"{out_dir}/{prefix}.EEs.tax.tsv")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.cl.tax",
              f"{out_dir}/{prefix}.EEs.cleaned.tax.tsv")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.fmt.fasta",
              f"{out_dir}/{prefix}.EEs.flanks.fa")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left.fasta",
              f"{out_dir}/{prefix}.EEs.L-flank.fa")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.left.fasta.tblastn.filtred",
              f"{out_dir}/{prefix}.EEs.L-flank.blast.tsv")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right.fasta",
              f"{out_dir}/{prefix}.EEs.R-flank.fa")
    os.rename(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.tax.bed.merge.fmt.fa.bed.flank.right.fasta.tblastn.filtred",
              f"{out_dir}/{prefix}.EEs.R-flank.blast.tsv")

    print("\n|"+"-"*45+"SUMMARY RESULTS"+"-"*45+"|\n")
    print("\n|"+"-"*45+"SUMMARY RESULTS"+"-"*45+"|\n", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.fa ----------------------------- Fasta file with Endogenous Elements nucleotide sequences.")
    print(f"{out_dir}/{prefix}.EEs.tax.tsv ------------------------ TSV file with Endogenous Elements taxonomy.")
    print(f"{out_dir}/{prefix}.EEs.flanks.fa ---------------------- Fasta file with Endogenous Elements plus {flank_region}nt in each flanking regions.")
    print(f"{out_dir}/{prefix}.EEs.L-flank.fa --------------------- Fasta file with Endogenous Elements plus {flank_region}nt upstream flanking region.")
    print(f"{out_dir}/{prefix}.EEs.R-flank.fa --------------------- Fasta file with Endogenous Elements plus {flank_region}nt downstream flanking region.")
    print(f"{out_dir}/{prefix}.EEs.L-flank.blast.tsv -------------- TSV file with filtred blast results of upstream flanking regions.")
    print(f"{out_dir}/{prefix}.EEs.R-flank.blast.tsv -------------- TSV file with filtred blast results of downstream flanking regions.")
    print(f"{out_dir}/{prefix}.EEs.cleaned.fa --------------------- Fasta file with Cleaned Endogenous Elements.")
    print(f"{out_dir}/{prefix}.EEs.cleaned.tax.tsv ---------------- TSV file with Cleaned Endogenous Elements.")

    print(f"{out_dir}/{prefix}.EEs.fa ----------------------------- Fasta file with Endogenous Elements nucleotide sequences.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.tax.tsv ------------------------ TSV file with Endogenous Elements taxonomy.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.flanks.fa ---------------------- Fasta file with Endogenous Elements plus {flank_region}nt in each flanking regions.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.L-flank.fa --------------------- Fasta file with Endogenous Elements plus {flank_region}nt upstream flanking region.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.R-flank.fa --------------------- Fasta file with Endogenous Elements plus {flank_region}nt downstream flanking region.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.L-flank.blast.tsv -------------- TSV file with filtred blast results of upstream flanking regions.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.R-flank.blast.tsv -------------- TSV file with filtred blast results of downstream flanking regions.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.cleaned.fa --------------------- Fasta file with Cleaned Endogenous Elements.", file=log_file)
    print(f"{out_dir}/{prefix}.EEs.cleaned.tax.tsv ---------------- TSV file with Cleaned Endogenous Elements.", file=log_file)

    if temp_remove == True:
        for tmp_file in glob.glob(f"{out_dir}/*rn*"):
            os.remove(tmp_file)
        print('\nTemporary files were removed.')
        print('\nTemporary files were removed.', file=log_file)
    else:
        if os.path.isdir(f'{out_dir}/tmp_files') == False:
            os.mkdir(f'{out_dir}/tmp_files')
        else:
            pass
        for tmp_file in glob.glob(f"{out_dir}/*rn*"):
            new_tmp_file = re.sub(r'.*/', '', tmp_file)
            os.rename(tmp_file, f"{out_dir}/tmp_files/{new_tmp_file}")
        print(
            f'\nTemporary files were moved to {out_dir}/tmp_files. Check the github documentation to access the description of each temporary file.')
        print(
            f'\nTemporary files were moved to {out_dir}/tmp_files. Check the github documentation to access the description of each temporary file.', file=log_file)
    print_info.print_finish()

    end_running_time = time.time()
    total_running_time = end_running_time - start_running_time
    print(f'TOTAL TIME = {total_running_time/60:.2f} MINUTES')
    print(f'TOTAL TIME = {total_running_time/60:.2f} MINUTES', file=log_file)
