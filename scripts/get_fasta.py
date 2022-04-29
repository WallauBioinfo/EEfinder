import subprocess
import shlex


def get_fasta(in_file, bed_file, out_file, log):
    """
    This function execute the bedtools getfasta.

    Keyword arguments:
    in_file - input_file, parsed with -in argument.
    bed_file - bed file, genereated along the pipeline
    out_file - output file 
    """

    get_fasta = f'bedtools getfasta -fi {in_file} -bed {bed_file} -fo {out_file}'
    get_fasta = shlex.split(get_fasta)
    cmd_get_fasta = subprocess.Popen(get_fasta)
    cmd_get_fasta.wait()
    print(f"DONE: Get fasta for {bed_file}!", file = log)
    return(print(f"DONE: Get fasta for {bed_file}!"))


class GetFasta:
    def __init__(self, in_file, bed_file, out_file, log):
        self.in_file = in_file
        self.bed_file = bed_file
        self.out_file = out_file
        self.log = log

    def run_get_fasta(self):
        get_fasta(self.in_file, self.bed_file, self.out_file, self.log)
