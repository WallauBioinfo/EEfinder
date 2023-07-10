import re


class InsertPrefix:
    """
    Insert prefix on input file fasta sequence headers

    Keyword arguments:
    input_file: input fasta file
    prefix: prefix string name
    outdir: output directory
    """

    def __init__(self, input_file: str, prefix: str, outdir: str) -> None:
        self.input_file = input_file
        self.prefix = prefix
        self.outdir = outdir

        self.insert_prefix()

    def insert_prefix(self) -> None:
        with open(self.input_file, "r") as input_file_reader, open(
            f"{self.outdir}/{self.prefix}.rn", "w+"
        ) as output_file:
            read_file = input_file_reader.read()
            for line in read_file:
                if ">" in line:
                    line = re.sub(r">", f">{self.prefix}/", line)
                output_file.write(line)
