import click

def print_paper_info():
    print("\n")
    print("Mensagem de Cabeçalho")
    print("\n")


def print_end_info():
    print("\n")
    print("|" + "-" * 48 + "FINISHING" + "-" * 48 + "|")
    print("\n")
    click.secho(f"Thank You for use EEfinder", fg="green")
    click.secho(f"Please Cite: Paper info", fg="green")
    print("\n")


class PaperInfo:
    def print_message(self):
        print_paper_info()

    def print_finish(self):
        print_end_info()
