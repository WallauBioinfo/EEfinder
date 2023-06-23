from Bio import SeqIO


def get_length(input_file, log):
    """
    This function creates a length file with the module SeqIO

    Keywords arguments:
    input_file = formated genome
    """
    with open(f"{input_file}.rn.fmt.lenght", "w") as output_length:
        length_list = []
        # writer_output_length = csv.writer(output_length,delimiter='\t')
        for seq_record in SeqIO.parse(input_file, "fasta"):
            output_length.write(f"{seq_record.id}\t{str(len(seq_record))}\n")
        print("DONE: Get Scaffold lengths!", file=log)
        return print("DONE: Get Scaffold lengths!")


class GetLength:
    def __init__(self, input_file, log):
        self.input_file = input_file
        self.log = log

    def run_get_length(self):
        get_length(self.input_file, self.log)
