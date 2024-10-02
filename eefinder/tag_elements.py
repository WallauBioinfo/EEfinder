import pandas as pd
import numpy as np


def _list_to_string(overlaped_elements: list) -> str:
    return ",".join(map(str, overlaped_elements))


class TagElements:
    """
    Create a collumn in the tax file with a tag sinalizing if the element is overlaped or unique
    
    Keyword arguments:
    tax_file: tax input file
    """
    def __init__(self, tax_file: str) -> object:
        self.tax_file = tax_file

        self.tag_elements()
        self.get_average_pident()

    def tag_elements(self) -> None:
        df = pd.read_csv(self.tax_file, sep="\t")
        df["Element-ID"] = df["Element-ID"].replace(".*/", "", regex=True)
        df["start"] = df["Element-ID"].replace("-.*", "", regex=True)
        df["start"] = df["start"].replace(".*:", "", regex=True)
        df["end"] = df["Element-ID"].replace(".*-", "", regex=True)
        df["contig"] = df["Element-ID"].replace(":.*", "", regex=True)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        matched_elements = []
        for i, row in df.iterrows():
            matched_ids = []
            filtered_df = df[df["contig"] == row["contig"]]
            for j, other_row in filtered_df.iterrows():
                if (
                    i != j
                    and other_row["start"] <= row["end"] + 100
                    and other_row["end"] >= row["start"] - 100
                    and other_row["Family"] != row["Family"]
                ):
                    matched_ids.append(other_row["Element-ID"])
            matched_elements.append(matched_ids)
        df["Overlaped_Element_ID"] = matched_elements

        df["Overlaped_Element_ID"] = df["Overlaped_Element_ID"].apply(_list_to_string)
        df.drop(columns=["start", "end", "contig"], inplace=True)

        for i, row in df.iterrows():
            if row["Overlaped_Element_ID"] == "":
                df.at[i, "tag"] = "unique"
            else:
                row["Overlaped_Element_ID"] != ""
                df.at[i, "tag"] = "overlaped"

        df.to_csv(self.tax_file, sep="\t", index=False, header=True)

    def get_average_pident(self) -> None:
        df = pd.read_csv(self.tax_file, sep="\t")
        def calculate_average(pidents):
            entries = pidents.split(" | ")
            pidents_values = [float(entry.split("|")[1]) for entry in entries if "|" in entry]
            return round(np.mean(pidents_values), 1) if pidents else np.nan
        df['Average_pident'] = df['Protein-IDs'].apply(calculate_average)
        df.to_csv(self.tax_file, sep="\t", index=False, header=True)