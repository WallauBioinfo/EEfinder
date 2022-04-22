import re


def insert_prefix(input_file, prefix, outdir):
    """
    This function insert the prefix in the fasta header

    Keyword arguments:
    input_file, parsed with -in argument.
    """
    with open(input_file, 'r') as input_file_reader, open(f"{outdir}/{prefix}.rn", 'w+') as output_file:
        read_file = input_file_reader.read()
        for line in read_file:
            if '>' in line:
                line = re.sub(r">", f">{prefix}/", line)
            output_file.write(line)


def concat_files(file1, file2, outdir, output):
    """
    This function concatenate two files.

    Keyword arguments:
    file1 - first file to be concatenated
    file2 - second file to be concatenated
    outdir - directory to store outputs
    output - output file name
    """
    with open(file1) as file1Read, open(file2) as file2Read:
        file1Read = file1Read.read()
        file2Read = file2Read.read()
        output_file = ''
        output_file += file1Read
        output_file += file2Read
        with open(f"{outdir}/{output}", 'w') as outputWriter:
            outputWriter.write(output_file)

        del file1Read, file2Read, output_file
    print(
        f"DONE: The concatenated filter database is stored at {outdir}{output}\n")


class SetPrefix():
    def __init__(self, in_file, prefix, outdir):
        self.in_file = in_file
        self.prefix = prefix
        self.outdir = outdir

    def run_insert_prefix(self):
        insert_prefix(self.in_file, self.prefix, self.outdir)


class ConcatFiles():
    def __init__(self, file1, file2, outdir, output):
        self.file1 = file1
        self.file2 = file2
        self.outdir = outdir
        self.output = output

    def run_concate_file(self):
        concat_files(self.file1, self.file2, self.outdir, self.output)
