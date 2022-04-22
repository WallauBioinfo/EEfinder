#!/usr/bin/python3
# -*- coding: utf-8 -*-

##################################>IMPORT-MODULES<##################################

import argparse
import re
from scripts.compare_results import CompareResults
from scripts.run_message import PaperInfo
from scripts.prepare_data import ConcatFiles, SetPrefix
from scripts.clean_data import RemoveShort, MaskClean
from scripts.make_database import MakeDB
from scripts.similarity_analysis import SimilaritySearch
from scripts.filter_table import FilterTable
from scripts.get_fasta import GetFasta
from scripts.get_taxonomy import GetTaxonomy
from scripts.format_bed import GetAnnotBed, RemoveAnnotation, ReformatBedName
from scripts.bed_merge import MergeBed
###############################>GENERAL-INFORMATIONS<###############################

__authors__ = {"Names": ["Filipe Dezordi", "Yago Dias"],
               "emails": ["zimmer.filipe@gmail.com", "yag.dias@gmail"],
               "gits": ["https://github.com/dezordi", "https://github.com/yagdias"]}
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Filipe Dezordi"
__status__ = "Prototype"

#####################################>ARGUMENTS<#####################################
parser = argparse.ArgumentParser(description='This tool predict regions of Endogenous Viral Elements.',
                                 usage='\n>Running with default parameters:\n   python EVEfinder.py -in <genome.fasta> -db <viral_prot.fasta> -mt <viral_prot.csv> -db1 <host_genes.fasta> -db2 <host_tranposons.fasta>. For more examples, see the documentation on Github.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "-in", "--input", help="Fasta file (nucleotides).")  # required=True
parser.add_argument("-od", "--outdir", help="Path and dir to store output")
parser.add_argument("-db", "--database",
                    help="Viral proteins database file.")  # required=True
parser.add_argument("-mt", "--dbmetadata",
                    help="Viral proteins csv file.")  # required=True
parser.add_argument("-db2", '--databaseHOST',
                    help="Host and heatshock proteins database file.")  # required=True
parser.add_argument("-db3", "--databaseTE",
                    help="Transposon protein database file.")  # required=True
parser.add_argument("-md", "--mode", help="Choose between Blast or the Diamond strategies (fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensisitve) to run analysis, default = blastx.",
                    default='blastx', choices=['blastx', 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'])
parser.add_argument(
    "-ln", "--length", help="Minimum length of contigs used for BLAST, default = 10,000.", type=int, default=10000)
parser.add_argument(
    "-fl", "--flank", help="Length of flanking regions of EVEs to be extracted, default = 10,000.", type=int, default=10000)
parser.add_argument(
    "-lm", "--limit", help="Limit of bases used to merge regions on bedtools merge, default = 1.", type=str, default=str(1))
parser.add_argument("-mp", "--mask_per",
                    help="Limit of lowercase letters in percentage to consider a putative EVE as a repetitive region, default = 50.", type=str, default=50)
parser.add_argument(
    "-p", "--threads", help="Threads for multi-thread analysis, default = 1.", type=int, default=1)
parser.add_argument("-rm", "--removetmp", help="Remove temporary files generated through analysis? default = True.",
                    default=False, choices=['True', 'False'])
parser.add_argument("-mk", "--makeblastdb", help="Make blast database?, default = True.",
                    default=True, choices=['True', 'False'])
parser.add_argument("-st", "--step", help="Select step of analysis, default = All.",
                    default='All', choices=['All', 'EVEfinder', 'Tax', 'TEs', 'Flank', 'Flank-blast'])
parser.add_argument("-pr", "--prefix", help="Write the prefix name for output files. We strongly recommend the use of genome name assembly, this prefix will be used to create the EVEname (The EVE name will be formated as PREFIX|CONTIG/SCAFFOLD:START-END) default = input file name.")

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
step = args.step
prefix = args.prefix

# Set Prefix and StorePath
if prefix == None:  # format prefix name, if the -pr argument was not used
    prefix = input_file
    prefix = re.sub("\..*", "", prefix)
    prefix = re.sub('.*/', '', prefix).rstrip('\n')
    prefix = re.sub('.*/', '', prefix).rstrip('\n')

# checar se n tiver um diretorio, criar um diretorio
if out_dir != None:
    if "/" in out_dir:
        out_dir = re.sub("/$", "", out_dir)
else:
    if '/' not in input_file:
        out_dir = '.'
    else:
        out_dir = re.sub("/.*", "", input_file).rstrip('\n')

if __name__ == '__main__':
    print_info = PaperInfo()
    print_info.print_message()
    # Prepate Data Step
    # Contenate filter files
    print("|"+"-"*45+"PREPARING DATA"+"-"*45+"|\n")
    concat_files = ConcatFiles(
        database_HOST, database_TE, out_dir, 'DatabaseFilter.fa')
    concat_files.run_concate_file()
    # Set Prefix
    set_prefix = SetPrefix(input_file, prefix, out_dir)
    set_prefix.run_insert_prefix()
    # Remove short sequences
    remove_seqs = RemoveShort(f"{out_dir}/{prefix}.rn", length_cutoff)
    remove_seqs.run_remove_short()
    # Make databases for EVEs screening and filters
    print("|"+"-"*42+"FORMATTING DATABASES"+"-"*42+"|\n")
    make_db_virus = MakeDB(mode, database_file, 'prot', threads)
    make_db_virus.run_make_db()
    make_db_filter = MakeDB(
        mode, f"{out_dir}/DatabaseFilter.fa", 'prot', threads)
    make_db_filter.run_make_db()
    # Run Similarity Screenings
    print("|"+"-"*40+"RUNNING SIMILARITY SEARCH"+"-"*39+"|\n")
    vir_similarity_step = SimilaritySearch(
        f"{out_dir}/{prefix}.rn.fmt", database_file, threads, mode)
    vir_similarity_step.run_similarity_search()
    # Filtering results
    vir_filter_table = FilterTable(f"{out_dir}/{prefix}.rn.fmt.blastx", "VIR")
    vir_filter_table.run_filter()
    # Getting fastas
    getter_fasta = GetFasta(f"{out_dir}/{prefix}.rn.fmt", f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed",
                            f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta")
    getter_fasta.run_get_fasta()
    # Second Similarity Screen

    host_similarity_step = SimilaritySearch(
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta", f"{out_dir}/DatabaseFilter.fa", threads, mode)
    host_similarity_step.run_similarity_search()
    # Filtering results
    host_filter_table = FilterTable(
        f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx", "HOST")
    host_filter_table.run_filter()
    # Comparing results
    comparer = CompareResults(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred",
                              f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred")
    comparer.run_comparation()
    # Get info taxonomy for pEVEs
    get_info = GetTaxonomy(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.vir",database_table)
    get_info.run_get_taxonomy_info()
    # Get annotated bed file
    get_annot_bed = GetAnnotBed(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.vir.tax")
    get_annot_bed.run_get_annotated_bed()
    # Merge bed file
    merge_bed = MergeBed(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.vir.tax.bed",str(limit_merge))
    merge_bed.run_merge_bedfile()
    # Format bed names to get_fasta step and to create a taxonomy table with merged elements
    remove_annot_bed = RemoveAnnotation(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.vir.tax.bed.merge")
    remove_annot_bed.run_reformat_bed()
    reformat_bed_name = ReformatBedName(f"{out_dir}/{prefix}.rn.fmt.blastx.filtred.bed.fasta.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt")
    reformat_bed_name.run_reformat_name()