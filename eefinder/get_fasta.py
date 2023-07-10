import subprocess
import shlex


class GetFasta:
    def __init__(self, input_file, bed_file, out_file):
        self.input_file = input_file
        self.bed_file = bed_file
        self.out_file = out_file

        self.get_fasta()

    def get_fasta(self):
        """
        This function execute the bedtools getfasta.

        Keyword arguments:
        input_file - input_file, parsed with -in argument.
        bed_file - bed file, genereated along the pipeline
        out_file - output file
        """

        get_fasta = f"bedtools getfasta -fi {self.input_file} -bed {self.bed_file} -fo {self.out_file}"
        get_fasta = shlex.split(get_fasta)
        cmd_get_fasta = subprocess.Popen(get_fasta)
        cmd_get_fasta.wait()
