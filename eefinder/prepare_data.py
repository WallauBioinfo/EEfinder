import re


def insert_prefix(input_file, prefix, outdir):
    """
    This function insert the prefix in the fasta header

    Keyword arguments:
    input_file, parsed with -in argument.
    """
    with open(input_file, "r") as input_file_reader, open(
        f"{outdir}/{prefix}.rn", "w+"
    ) as output_file:
        read_file = input_file_reader.read()
        for line in read_file:
            if ">" in line:
                line = re.sub(r">", f">{prefix}/", line)
            output_file.write(line)


class SetPrefix:
    def __init__(self, in_file, prefix, outdir):
        self.in_file = in_file
        self.prefix = prefix
        self.outdir = outdir

    def run_insert_prefix(self):
        insert_prefix(self.in_file, self.prefix, self.outdir)
