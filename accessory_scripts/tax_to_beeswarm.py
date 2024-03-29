#!/usr/bin/python3
# -*- coding: utf-8 -*-
###############################>GENERAL-INFORMATIONS<###############################
"""
Build in Python 3.6.5+

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://dezordi.github.io/

Script repository:
https://github.com/dezordi/PEVEI

"""
###############################>LIBRARIES<###############################
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.style as style

###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(description="Create beeswarm plots from .tax files")
parser.add_argument(
    "-in", "--input", help=".Tax file generated by PEVEI.py", required=True
)
parser.add_argument(
    "-md",
    "--mode",
    help="Create plots with family or order taxonomy? default=Family",
    default="Family",
    choices=["Family", "Order"],
)
parser.add_argument(
    "-st",
    "--specifictax",
    help="Pass specific taxon to generate plots, e.g. Rhabdoviridae Flaviviridae",
    nargs="+",
)
parser.add_argument(
    "-ms",
    "--markersize",
    help="Pass the marker sizer length (dots on beeswarm plot)",
    default=6,
    type=int,
)
parser.add_argument(
    "-gt",
    "--groupbytax",
    help="Group EVEs families/orders by the number of elements, default = 5",
    default=5,
    type=int,
)
parser.add_argument(
    "-gp",
    "--groupbyprot",
    help="Group EVEs proteins categories by the number of proteins, default = 20 ",
    default=20,
    type=int,
)
args = parser.parse_args()
input_file = args.input
plot_mode = args.mode
specific_tax = args.specifictax
mark_size = args.markersize
group_bytax = args.groupbytax
group_byprot = args.groupbyprot
###############################>SNS-STYLE<###############################
sns.set()
style.use("seaborn-poster")
sns.set_palette("husl", 9)
###############################>DATAFRAME<###############################
df = pd.read_csv(input_file, sep="\t", header=0)
if specific_tax != None:
    df = df.loc[df["Family"].isin(specific_tax)]

family_count = df["Family"].value_counts()
family_low = family_count[family_count <= group_bytax]
family_low_list = family_low.index
for i in family_low_list:
    df["Family"] = df["Family"].replace({i: "Others"})
# grouping proteins terms with less than 5 elements to 'Other proteins'
protein_count = df["ProteinProducts"].value_counts()
protein_low = protein_count[protein_count <= group_byprot]
protein_low_list = protein_low.index
for i in protein_low_list:
    df["ProteinProducts"] = df["ProteinProducts"].replace({i: "Other proteins"})
# get number of elements by strand setnse
sense_count = df["Sense"].value_counts()
pos_count = str(sense_count["pos"])
neg_count = str(sense_count["neg"])
###############################>GENERAL-BEESWARNPLOT<###############################
beeswarn_plot = sns.swarmplot(
    x="Length", y="Super-Kingdom", hue=plot_mode, data=df, size=mark_size
)
# sns.despine(fig=None, top=True, right=True, left=True, bottom=False, offset=None, trim=False)
beeswarn_plot.set_ylabel("")
beeswarn_plot.set_yticks([])
beeswarn_plot.set_xlabel("EVE length (pb)")
beeswarn_plot.legend(
    loc="upper right",
    ncol=2,
    fancybox=True,
    prop={"size": 12},
    title="EVEs " + plot_mode,
    title_fontsize=14,
)
beeswarn_plot.set_xticks([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
beeswarn_plot.spines["bottom"].set_color("0")
beeswarn_plot.spines["top"].set_color("0")
beeswarn_plot.spines["right"].set_color("0")
beeswarn_plot.spines["left"].set_color("0")
plt.savefig(input_file + ".beeswarn_plot.pdf", dpi=300)
plt.clf()
###############################>SENSE-BEESWARNPLOT<###############################
beeswarn_plot = sns.swarmplot(
    x="Sense", y="Length", hue=plot_mode, data=df, size=mark_size
)
# sns.despine(fig=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
beeswarn_plot.set_xlabel("Strand Sense")
beeswarn_plot.set_ylabel("EVE length (pb)")
beeswarn_plot.set_xticklabels(
    ["Negative (n = " + neg_count + ")", "Positive (n = " + pos_count + ")"]
)
beeswarn_plot.legend(
    loc="upper center",
    ncol=2,
    fancybox=True,
    prop={"size": 12},
    title="EVEs " + plot_mode,
    title_fontsize=14,
)
beeswarn_plot.set_yticks([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
beeswarn_plot.spines["bottom"].set_color("0")
beeswarn_plot.spines["top"].set_color("0")
beeswarn_plot.spines["right"].set_color("0")
beeswarn_plot.spines["left"].set_color("0")
plt.savefig(input_file + ".beeswarn_plot_strand.pdf", dpi=300, bbox_inches="tight")
plt.clf()
###############################>PROTEIN-BEESWARNPLOT<###############################
beeswarn_plot = sns.swarmplot(
    x="Length", y="ProteinProducts", hue=plot_mode, data=df, size=mark_size
)
# sns.despine(fig=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
beeswarn_plot.set_ylabel("Protein term")
beeswarn_plot.set_xlabel("EVE length (pb)")
beeswarn_plot.legend(
    loc="upper right",
    ncol=2,
    fancybox=True,
    prop={"size": 12},
    title="EVEs " + plot_mode,
    title_fontsize=14,
)
beeswarn_plot.set_xticks([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
beeswarn_plot.spines["bottom"].set_color("0")
beeswarn_plot.spines["top"].set_color("0")
beeswarn_plot.spines["right"].set_color("0")
beeswarn_plot.spines["left"].set_color("0")
plt.savefig(input_file + ".beeswarn_plot_product.pdf", dpi=300, bbox_inches="tight")
