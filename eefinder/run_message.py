import click


def print_paper_info(version: str):
    print("\n")
    print("  ______ ______  __ _           _           ")
    print(" |  ____|  ____|/ _(_)         | |          ")
    print(" | |__  | |__  | |_ _ _ __   __| | ___ _ __ ")
    print(" |  __| |  __| |  _| | '_ \ / _` |/ _ \ '__|")
    print(" | |____| |____| | | | | | | (_| |  __/ |   ")
    print(" |______|______|_| |_|_| |_|\__,_|\___|_|   ")
    print(f"                            version {version}")
    print("\n")


def print_end_info():
    click.secho("Thank You for use EEfinder!", fg="green")
    click.secho("Please Cite: Dias, Y. J. M., Dezordi, F. Z., & Wallau, G. L. (2024). "
                +"EEFinder: A general-purpose tool for identification of bacterial and viral endogenized elements in eukaryotic genomes. "
                +"Computational and Structural Biotechnology Journal. https://doi.org/10.1016/j.csbj.2024.10.012", fg="green")
    print("\n")


class PaperInfo:
    def print_start(self, version: str):
        print_paper_info(version)

    def print_finish(self):
        print_end_info()
