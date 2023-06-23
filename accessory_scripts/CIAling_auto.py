#!/usr/bin/python3
# -*- coding: utf-8 -*-
###############################>GENERAL-INFORMATIONS<###############################
"""
Build in Python 3.6

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://github.com/dezordi

Dependencies:
CD-HIT version 4.7
"""
###############################>LIBRARIES<###############################

import argparse, csv, os, subprocess, shlex, sys, time, re
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    description="This scripts receives a fasta file and execute a cd-hit-est analysis and formats the cd-hit output",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument("-in", "--input", help="Fasta file", required=True)
parser.add_argument(
    "-p", "--threads", help="Threads for Mafft analysis", default=1, type=str
)
# Storing argument on variables
args = parser.parse_args()
input_file = args.input
var_threads = args.threads
# Creating directory for brute
dir_name = re.sub(r"\.", "_", input_file)
dir_name = re.sub(r"fasta", "CIAlign_results", dir_name)
create_dir = subprocess.call([f"mkdir {dir_name}"], shell=True)
# Get average distance
def get_avg(tsv_file):
    with open(tsv_file, "r") as matrix_file:
        matrix_df = pd.read_csv(matrix_file, sep="\t", index_col=0)
        matrix_df.replace(1, np.nan, inplace=True)
        mean_df = matrix_df.mean()
        mean_df_list = mean_df.tolist()
        matrix_file.close()
    global average_similarity
    average_similarity = "{:.3f}".format(sum(mean_df_list) / len(mean_df_list))


# Align sequences
with open(input_file + ".algn", "w") as mafft_out:
    mafft_cmd = (
        f"mafft --thread "
        + str(var_threads)
        + " --reorder --adjustdirection --bl 45 --auto "
        + input_file
    )
    mafft_cmd = shlex.split(mafft_cmd)
    mafft_cmd_process = subprocess.Popen(mafft_cmd, stdout=mafft_out)
    mafft_cmd_process.wait()
# Make a 1s distance matrix
CIAling_matrix = subprocess.call(
    [
        f"CIAlign --infile {input_file}.algn --outfile_stem {input_file}.matrix1 --make_similarity_matrix_input"
    ],
    shell=True,
)
# Get 1st distan treshold
get_avg(input_file + ".matrix1_input_similarity.tsv")
print(f"first distance: {average_similarity}")
# Removing sequences using the 1st distance matrix treshold
CIAling_matrix = subprocess.call(
    [
        f"CIAlign --infile {input_file}.algn --outfile_stem {input_file}.01divergent --remove_divergent --remove_divergent_minperc {average_similarity} --plot_markup"
    ],
    shell=True,
)
# Remove Inserions
CIAling_matrix = subprocess.call(
    [
        f"CIAlign --infile {input_file}.01divergent_cleaned.fasta --outfile_stem {input_file}.02insertions --remove_insertions --plot_markup"
    ],
    shell=True,
)
# Crop ends
CIAling_matrix = subprocess.call(
    [
        f"CIAlign --infile {input_file}.02insertions_cleaned.fasta --outfile_stem {input_file}.03ends --crop_ends --plot_markup"
    ],
    shell=True,
)
# Make a 2nd distance matrix
CIAling_matrix = subprocess.call(
    [
        f"CIAlign --infile {input_file}.03ends_cleaned.fasta --outfile_stem {input_file}.matrix2 --make_similarity_matrix_input"
    ],
    shell=True,
)
# Get 2nd distan treshold
get_avg(input_file + ".matrix2_input_similarity.tsv")
print(f"second distance: {average_similarity}")
# Removing sequences using the 2nd distance matrix treshold
CIAling_matrix = subprocess.call(
    [
        f"CIAlign --infile {input_file}.03ends_cleaned.fasta --outfile_stem {input_file}.04divergent2 --remove_divergent --remove_divergent_minperc {average_similarity} --plot_markup"
    ],
    shell=True,
)
# Organising results
mv_output = subprocess.call(
    [f"mv {input_file}.04divergent2_cleaned.fasta {input_file}.algn.edited"], shell=True
)
move_results = subprocess.call(
    [f"mv *divergent* *matrix* *insertions* *ends* {dir_name}"], shell=True
)
