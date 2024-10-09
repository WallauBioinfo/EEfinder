import pandas as pd
import csv
import os
import re
import glob
import shutil


class FilterTable:
    """
    Receives a blastx result and filter based on query ID and ranges of qstart and qend.

    Keyword arguments:
    blast_result: input blastx result
    rangejunction: range for filter redundant hits
    tag: HOST or EE, tells which blastx is
    out_dir: output directory, parsed by -od
    """

    def __init__(
        self, blast_result: str, rangejunction: int, tag: str, out_dir: str
    ) -> object:
        self.blast_result = blast_result
        self.rangejunction = rangejunction
        self.tag = tag
        self.out_dir = out_dir

        self.filter_blast()

    def filter_blast(self) -> None:
        header_outfmt6 = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ]  # creates a blast header output in format = 6
        df = pd.read_csv(
            self.blast_result, sep="\t", header=None, names=header_outfmt6
        ).sort_values(by="bitscore", ascending=False)
        df["sense"] = ""
        df["bed_name"] = ""
        df["tag"] = ""
        df["new_qstart"] = df["qstart"]
        df["new_qend"] = df["qend"]
        df.to_csv(self.blast_result + ".csv", sep="\t")
        chunks = df = pd.read_csv(
            f"{self.blast_result}.csv", sep="\t", chunksize=200000
        )
        count = 0
        tmp_path = f"{self.out_dir}/tmp/"
        if os.path.exists(tmp_path) == False:
            os.mkdir(tmp_path)
        for df in chunks:
            df["sense"] = df["sense"].astype(object)
            df.loc[
                df["qstart"].astype(int) > df["qend"].values.astype(int), "sense"
            ] = "neg"
            df.loc[
                df["qend"].values.astype(int) > df["qstart"].astype(int), "sense"
            ] = "pos"
            df.loc[df["sense"] == "neg", "new_qstart"] = df["qend"]
            df.loc[df["sense"] == "neg", "new_qend"] = df["qstart"]
            df.loc[df["sense"] == "neg", "qstart"] = df["new_qstart"]
            df.loc[df["sense"] == "neg", "qend"] = df["new_qend"]
            df.drop(columns=["new_qstart", "new_qend"], inplace=True)
            if self.tag == "EE":
                df["tag"] = "EE"
                df["bed_name"] = df.apply(
                    lambda x: "%s:%s-%s" % (x["qseqid"], x["qstart"], x["qend"]), axis=1
                )
            else:
                df["tag"] = "HOST"
                df["bed_name"] = df["qseqid"]
            pd.options.display.float_format = "{:,.2f}".format
            df["evalue"] = pd.to_numeric(df["evalue"], downcast="float")
            df = df[df.length >= 33]
            header = [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "sense",
                "bed_name",
                "tag",
            ]
            df = df[header]
            with open(f"{tmp_path}chunk.{count}.tsv", "w") as chunk_writer:
                df.to_csv(chunk_writer, sep="\t", index=False)
            count += 1
        all_chunks = glob.glob(f"{tmp_path}/*.tsv")
        final_filtred_file = pd.DataFrame()
        chunks_list = []
        for chunk in all_chunks:
            df = pd.read_csv(chunk, sep="\t")
            chunks_list.append(df)
        final_filtred_file = pd.concat(chunks_list, ignore_index=True)
        final_filtred_file["qstart_rng"] = final_filtred_file.qstart.floordiv(
            self.rangejunction
        )
        final_filtred_file["qend_rng"] = final_filtred_file.qend.floordiv(
            self.rangejunction
        )
        final_filtred_file = (
            final_filtred_file.drop_duplicates(subset=["qseqid", "qstart_rng", "sense"])
            .drop_duplicates(subset=["qseqid", "qstart_rng", "sense"])
            .sort_values(by=["qseqid"])
        )
        final_filtred_file.to_csv(
            f"{self.blast_result}.filtred", sep="\t", index=False, columns=header
        )
        if self.tag == "EE":
            final_filtred_file.to_csv(
                f"{self.blast_result}.filtred.bed",
                header=False,
                sep="\t",
                index=False,
                columns=["qseqid", "qstart", "qend"],
            )
        shutil.rmtree(tmp_path, ignore_errors=True)
