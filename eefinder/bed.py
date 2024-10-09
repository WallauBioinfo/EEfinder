import pandas as pd
import numpy as np
import shlex
import subprocess
import re


class GetFasta:
    """
    This function execute the bedtools getfasta.

    Keyword arguments:
    input_file: input_file, parsed with -in argument.
    bed_file: bed file, genereated along the pipeline
    out_file: output file
    """

    def __init__(self, input_file: str, bed_file: str, out_file: str) -> object:
        self.input_file = input_file
        self.bed_file = bed_file
        self.out_file = out_file

        self.get_fasta()

    def get_fasta(self) -> None:
        get_fasta = f"bedtools getfasta -fi {self.input_file} -bed {self.bed_file} -fo {self.out_file}"
        get_fasta = shlex.split(get_fasta)
        cmd_get_fasta = subprocess.Popen(
            get_fasta, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        cmd_get_fasta.wait()


class GetAnnotBed:
    """
    Create a bed file that will be used to merge truncated EVEs of the same
    family in the same sense based on a limite length treshold.

    Keyword arguments:
    blast_tax_info: csv file generated in the get_taxonomy_info function on get_taxonomy.py
    merge_level: genus or family, choose which level going to merge nearby elements
    """

    def __init__(self, blast_tax_info: str, merge_level: str) -> object:
        self.blast_tax_info = blast_tax_info
        self.merge_level = merge_level

        self.get_annotated_bed()

    def get_annotated_bed(self) -> None:
        df_blast_tax_info = pd.read_csv(self.blast_tax_info, sep=",")
        df_blast_tax_info["qseqid"] = df_blast_tax_info["qseqid"].str.replace(
            r"\:.*", "", regex=True
        )
        df_blast_tax_info["sseqid"] = (
            df_blast_tax_info["sseqid"]
            + "|"
            + df_blast_tax_info["sense"]
            + "|"
            + df_blast_tax_info["pident"].astype(str)
        )
        df_blast_tax_info["Family"] = df_blast_tax_info["Family"].fillna("Unknown")
        df_blast_tax_info["Genus"] = df_blast_tax_info["Genus"].fillna("Unknown")

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
    """
    Remove the annotated information generate into the get_annotated_bed function.

    Keyword arguments:
    bed_annotated_merged_file: tsv file generated in the merge_bedfile function on bed_merge.py
    """

    def __init__(self, bed_annotated_merged_file: str) -> object:
        self.bed_annotated_merged_file = bed_annotated_merged_file

        self.reformat_bed()

    def reformat_bed(self) -> None:
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
    """
    Execute the bedtools merge.

    Keyword arguments:
    bed_annotated_file: annotated bed file created at get_annotated_bed function
    limit_merge: Limit of bases to merge regions, parsed with -lm argument
    """

    def __init__(self, bed_annotated_file: str, limit_merge: int) -> object:
        self.bed_annotated_file = bed_annotated_file
        self.limit_merge = limit_merge

        self.merge_bed()

    def merge_bed(self) -> None:
        bed_merge_output = open(f"{self.bed_annotated_file}.merge", "w")
        bed_merge_cmd = f'bedtools merge -d {int(self.limit_merge)} -i {self.bed_annotated_file} -c 4 -o collapse -delim " AND "'
        bed_merge_cmd = shlex.split(bed_merge_cmd)
        bed_merge_process = subprocess.Popen(bed_merge_cmd, stdout=bed_merge_output)
        bed_merge_process.wait()


class BedFlank:
    """
    Extract flanking regions of EEs using bedtools slop.

    Keyword arguments:
    input_file: bed file generated by get_bed function
    lenght_file: lenght file produced by get_length function
    flank_region: desired lenght regions for extraction, parsed from
    """

    def __init__(self, input_file: str, length_file: str, flank_region: int) -> object:
        self.input_file = input_file
        self.length_file = length_file
        self.flank_region = flank_region

        self.bedtools_flank()

    def bedtools_flank(self) -> None:
        with open(f"{self.input_file}.flank", "w") as flank_out:
            bed_flank_cmd = f"bedtools slop -i {self.input_file} -g {self.length_file} -b {str(self.flank_region)}"
            bed_flank_cmd = shlex.split(bed_flank_cmd)
            bed_flank_process = subprocess.Popen(bed_flank_cmd, stdout=flank_out)
            bed_flank_process.wait()


class GetBed:
    """
    Create a bed file from fasta file using replace logic.

    Keyword arguments:
    input_file: fasta file for desired bed file
    """

    def __init__(self, input_file: str) -> object:
        self.input_file = input_file

        self.get_bed()

    def get_bed(self) -> None:
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
