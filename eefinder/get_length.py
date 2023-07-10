from Bio import SeqIO


class GetLength:
    def __init__(self, input_file):
        self.input_file = input_file

        self.get_length()

    def get_length(self):
        """
        This function creates a length file with the module SeqIO

        Keywords arguments:
        input_file = formated genome
        """
        with open(f"{self.input_file}.rn.fmt.lenght", "w") as output_length:
            length_list = []
            for seq_record in SeqIO.parse(self.input_file, "fasta"):
                output_length.write(f"{seq_record.id}\t{str(len(seq_record))}\n")
