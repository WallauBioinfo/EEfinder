import pandas as pd
import numpy as np
import shlex
import subprocess
import re


class GetAnnotBed:
    def __init__(self, blast_tax_info, merge_level):
        self.blast_tax_info = blast_tax_info
        self.merge_level = merge_level

        self.get_annotated_bed()

    def get_annotated_bed(self):
        """
        This function create a bed file that will be used to merge truncated EVEs of the same family in the same sense based on a limite length treshold.

        Keyword arguments:
        blast_tax_info = csv file generated in the get_taxonomy_info function on get_taxonomy.py
        """

        df_blast_tax_info = pd.read_csv(self.blast_tax_info, sep=",")
        df_blast_tax_info["qseqid"] = df_blast_tax_info["qseqid"].str.replace(
            r"\:.*", "", regex=True
        )
        df_blast_tax_info["sseqid"] = (
            df_blast_tax_info["sseqid"] + "|" + df_blast_tax_info["sense"]
        )
        df_blast_tax_info["Family"].fillna("Unknown", inplace=True)
        df_blast_tax_info["Genus"].fillna("Unknown", inplace=True)

        if self.merge_level == "genus":
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
        bed_blast_info.to_csv(
            f"{self.blast_tax_info}.bed", index=False, header=False, sep="\t"
        )


class RemoveAnnotation:
    def __init__(self, bed_annotated_merged_file):
        self.bed_annotated_merged_file = bed_annotated_merged_file

        self.reformat_bed()

    def reformat_bed(self):
        """
        This function remove the annotated information generate into the get_annotated_bed function

        Keyword arguments:
        bed_annotated_merged_file = tsv file generated in the merge_bedfile function on bed_merge.py
        """

        df_merge_file = pd.read_csv(
            self.bed_annotated_merged_file, sep="\t", header=None
        )
        df_merge_file.iloc[:, 0] = df_merge_file.iloc[:, 0].str.replace(
            "\|.*", "", regex=True
        )
        df_merge_file.to_csv(
            f"{self.bed_annotated_merged_file}.fmt", index=False, header=False, sep="\t"
        )


class MergeBed:
    def __init__(self, bed_annotated_file, limit_merge):
        self.bed_annotated_file = bed_annotated_file
        self.limit_merge = limit_merge

        self.merge_bed()

    def merge_bed(self):
        bed_merge_output = open(f"{self.bed_annotated_file}.merge", "w")
        bed_merge_cmd = f'bedtools merge -d {self.limit_merge} -i {self.bed_annotated_file} -c 4 -o collapse -delim " AND "'
        bed_merge_cmd = shlex.split(bed_merge_cmd)
        bed_merge_process = subprocess.Popen(bed_merge_cmd, stdout=bed_merge_output)
        bed_merge_process.wait()


class BedFlank:
    def __init__(self, input_file, length_file, flank_region):
        self.input_file = input_file
        self.length_file = length_file
        self.flank_region = flank_region

        self.bedtools_flank()

    def bedtools_flank(self):
        with open(f"{self.input_file}.flank", "w") as flank_out:
            bed_flank_cmd = f"bedtools slop -i {self.input_file} -g {self.length_file} -b {str(self.flank_region)}"
            bed_flank_cmd = shlex.split(bed_flank_cmd)
            bed_flank_process = subprocess.Popen(bed_flank_cmd, stdout=flank_out)
            bed_flank_process.wait()


class GetBed:
    def __init__(self, input_file):
        self.input_file = input_file

        self.get_bed()

    def get_bed(self):
        with open(f"{self.input_file}", "r") as repeat_eves, open(
            f"{self.input_file}.bed", "w"
        ) as repeat_eves_bed_out:
            repeat_eves_lines = repeat_eves.readlines()
            for line in repeat_eves_lines:
                if ">" in line:
                    line_name = line.replace(">", "")
                    line_name = re.sub(":.*", "", line_name).rstrip("\n")
                    line_start = re.sub(".*:", "", line)
                    line_start = re.sub("-.*", "", line_start).rstrip("\n")
                    line_end = re.sub(".*-", "", line).rstrip("\n")
                    repeat_eves_bed_out.write(
                        f"{line_name}\t{line_start}\t{line_end}\n"
                    )
