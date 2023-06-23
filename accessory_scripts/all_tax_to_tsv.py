#!/usr/bin/python3
# -*- coding: utf-8 -*-
###############################>GENERAL-INFORMATIONS<###############################
"""
Build in Python 3.6

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://dezordi.github.io/

Script repository:
https://github.com/dezordi/PEVEI

Test files:
https://github.com/dezordi/PEVEI/test_files/all_tax_to_tsv/
"""
###############################>LIBRARIES<###############################

import pandas as pd
import numpy as np
import argparse, csv, os, subprocess, shlex, sys, time, re

###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(
    description="This scripts receives a list file with taxonomy results of PEVEI, and returns a tsv file with the taxonomy organized by assembly",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "-in", "--input", help="List with taxonomy file names. e.g.: tax.lst", required=True
)
# Storing argument on variables
args = parser.parse_args()
input_file = args.input
###############################>EXECUTION<###############################
# Creating lists to store genome assemblies and viral taxonomy
assembly_list = []
taxonomy_list = []
with open(input_file, "r") as tax_list:
    for tax_file in tax_list:
        tax_file = tax_file.rstrip("\n")
        with open(tax_file, "r") as tax_file_reader:
            tax_file_reader_csv = csv.reader(tax_file_reader, delimiter="\t")
            for line in tax_file_reader_csv:
                assembly = re.sub(r"/.*", "", line[0]).rstrip("\n")
                if assembly not in assembly_list:
                    assembly_list.append(assembly)
                taxonomy = "-".join([line[5], line[6]]).rstrip("\n")
                if taxonomy not in taxonomy_list:
                    taxonomy_list.append(taxonomy)
taxonomy_list.remove("Order-Family")
assembly_list.remove("Element-ID")
# Create a dataframe with viral taxonomy as index and assemblies code as columns
df = pd.DataFrame(index=taxonomy_list, columns=assembly_list)
with open(input_file, "r") as tax_list:
    for tax_file in tax_list:
        tax_file = tax_file.rstrip("\n")
        with open(tax_file, "r") as tax_file_reader:
            df_tax = pd.read_csv(tax_file_reader, sep="\t")
            taxonomy_formated = df_tax[["Order", "Family"]].agg("-".join, axis=1)
            df_tax["Taxonomy"] = taxonomy_formated
            df_tax = df_tax.replace(to_replace="/.*", value="", regex=True)
            for column in df.columns[0:]:
                if column in df_tax["Element-ID"].values:
                    for i, row in df.iterrows():
                        column_count = df_tax.Taxonomy.str.count(i).sum()
                        df[column].loc[i] = column_count
df.to_csv(input_file + ".tsv", sep="\t", na_rep="0")  # Export dataframe as tsv file
