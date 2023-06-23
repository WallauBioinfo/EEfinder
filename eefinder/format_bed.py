import pandas as pd
import numpy as np


def get_annotated_bed(blast_tax_info, merge_level, log):
    """
    This function create a bed file that will be used to merge truncated EVEs of the same family in the same sense based on a limite length treshold.

    Keyword arguments:
    blast_tax_info = csv file generated in the get_taxonomy_info function on get_taxonomy.py
    """

    df_blast_tax_info = pd.read_csv(blast_tax_info, sep=",")
    df_blast_tax_info["qseqid"] = df_blast_tax_info["qseqid"].str.replace(
        r"\:.*", "", regex=True
    )
    df_blast_tax_info["sseqid"] = (
        df_blast_tax_info["sseqid"] + "|" + df_blast_tax_info["sense"]
    )
    df_blast_tax_info["Family"].fillna("Unknown", inplace=True)
    df_blast_tax_info["Genus"].fillna("Unknown", inplace=True)
    """
    
    """
    if merge_level == "genus":
        df_blast_tax_info["formated_name"] = np.where(
            df_blast_tax_info["Genus"] != "Unknown",
            df_blast_tax_info["qseqid"]
            + "|"
            + df_blast_tax_info["Family"]
            + "|"
            + df_blast_tax_info["Genus"]
            + "|"
            + df_blast_tax_info["sense"],
            df_blast_tax_info["qseqid"]
            + "|"
            + df_blast_tax_info["sseqid"]
            + "|"
            + df_blast_tax_info["Genus"],
        )
    else:
        df_blast_tax_info["formated_name"] = np.where(
            df_blast_tax_info["Family"] != "Unknown",
            df_blast_tax_info["qseqid"]
            + "|"
            + df_blast_tax_info["Family"]
            + "|"
            + df_blast_tax_info["sense"],
            df_blast_tax_info["qseqid"]
            + "|"
            + df_blast_tax_info["sseqid"]
            + "|"
            + df_blast_tax_info["Family"],
        )
    bed_blast_info = df_blast_tax_info[
        ["formated_name", "qstart", "qend", "sseqid"]
    ].copy()
    bed_blast_info = bed_blast_info.sort_values(
        ["formated_name", "qstart"], ascending=(True, True)
    )
    bed_blast_info.to_csv(f"{blast_tax_info}.bed", index=False, header=False, sep="\t")
    print(f"DONE: Create annotated bed file.", file=log)
    return print(f"DONE: Create annotated bed file.")


def reformat_bed(bed_annotated_merged_file, log):
    """
    This function remove the annotated information generate into the get_annotated_bed function

    Keyword arguments:
    bed_annotated_merged_file = tsv file generated in the merge_bedfile function on bed_merge.py
    """

    df_merge_file = pd.read_csv(bed_annotated_merged_file, sep="\t", header=None)
    df_merge_file.iloc[:, 0] = df_merge_file.iloc[:, 0].str.replace(
        "\|.*", "", regex=True
    )
    df_merge_file.to_csv(
        f"{bed_annotated_merged_file}.fmt", index=False, header=False, sep="\t"
    )
    print(f"DONE: Remove annotation of bed file sequences name.", file=log)
    return print(f"DONE: Remove annotation of bed file sequences name.")


class GetAnnotBed:
    def __init__(self, blast_tax_info, merge_level, log):
        self.blast_tax_info = blast_tax_info
        self.merge_level = merge_level
        self.log = log

    def run_get_annotated_bed(self):
        get_annotated_bed(self.blast_tax_info, self.merge_level, self.log)


class RemoveAnnotation:
    def __init__(self, bed_annotated_merged_file, log):
        self.bed_annotated_merged_file = bed_annotated_merged_file
        self.log = log

    def run_reformat_bed(self):
        reformat_bed(self.bed_annotated_merged_file, self.log)
