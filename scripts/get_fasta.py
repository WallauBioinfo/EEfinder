import subprocess
import shlex


def get_fasta(in_file, bed_file, out_file):
    """
    This function execute the bedtools getfasta.

    Keyword arguments:
    in_file - input_file, parsed with -in argument.
    bed_file - bed file, genereated along the pipeline
    out_file - output file 
    """

    get_fasta = 'bedtools getfasta -fi '+in_file+' -bed '+bed_file+' -fo '+out_file
    get_fasta = shlex.split(get_fasta)
    cmd_get_fasta = subprocess.Popen(get_fasta)
    cmd_get_fasta.wait()
    print("DONE GETTING FASTA")


class GetFasta:
    def __init__(self, in_file, bed_file, out_file):
        self.in_file = in_file
        self.bed_file = bed_file
        self.out_file = out_file

    def run_get_fasta(self):
        get_fasta(self.in_file, self.bed_file, self.out_file)
