import pandas as pd


class CompareResults:
    """
    This function compares 2 blast results, for queries with same ID, only the one with the major bitscore is keept. In a final step
    Only queries with tag EE are maintained

    Keyword arguments:
    vir_result: filtred blast against ee database
    host_result: filtred blast against filter database
    """
    def __init__(self, vir_result: str, host_result: str) -> object:
        self.vir_result = vir_result
        self.host_result = host_result

        self.compare_results()

    def compare_results(self) -> None:
        df_vir = pd.read_csv(self.vir_result, sep="\t")
        df_vir["qseqid"] = df_vir["bed_name"]
        df_host = pd.read_csv(self.host_result, sep="\t")
        df_hybrid = pd.concat([df_vir, df_host], ignore_index=True)
        df_hybrid = df_hybrid.sort_values(by=["qseqid", "bitscore"], ascending=False)
        df_hybrid.to_csv(self.host_result + ".concat", sep="\t", index=False)
        df_nr = df_hybrid.drop_duplicates(subset=["qseqid"])
        df_nr.to_csv(self.host_result + ".concat.nr", sep="\t", index=False)
        df_nr_vir = df_nr[df_nr.tag == "EE"]
        df_nr_vir.to_csv(self.host_result + ".concat.nr", sep="\t", index=False)
