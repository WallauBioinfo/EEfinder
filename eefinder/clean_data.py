from Bio import SeqIO


class RemoveShortSequences:
    """
    Remove sequences bellow the cutoff threshold.

    Keyword arguments:
    input_file: input fasta file
    cutoff: cutoff length, parsed by -ln
    """

    def __init__(self, input_file: str, cutoff: str) -> object:
        self.input_file = input_file
        self.cutoff = cutoff

        self.cut_seq()

    def cut_seq(self) -> None:
        new_sequences = []
        input_handle = open(self.input_file, "r")
        output_handle = open(self.input_file + ".fmt", "w")
        for record in SeqIO.parse(input_handle, "fasta"):
            if len(record.seq) >= int(self.cutoff):
                new_sequences.append(record)
        SeqIO.write(new_sequences, output_handle, "fasta")


class MaskClean:
    """
    Remove sequences of EE on regions with a certain % of soft masked bases.

    Keyword arguments:
    input_file: fasta file, with putative EEs
    m_per: treshold masked percentage value, parsed with -mp argument
    """

    def __init__(self, input_file: str, m_per:str) -> object:
        self.input_file = input_file
        self.m_per = m_per

        self.mask_clean()

    def mask_clean(self) -> None:
        sequences = {}
        for seq_record in SeqIO.parse(self.input_file, "fasta"):
            sequence = str(seq_record.seq)
            sequence_id = str(seq_record.id)
            if (
                float(
                    sequence.count("a")
                    + sequence.count("t")
                    + sequence.count("c")
                    + sequence.count("g")
                    + sequence.count("n")
                    + sequence.count("N")
                )
                / float(len(sequence))
            ) * 100 <= float(self.m_per):
                if sequence_id not in sequences:
                    sequences[sequence_id] = sequence
        with open(self.input_file + ".cl", "w+") as output_file:
            for sequence_id, sequence in sequences.items():
                output_file.write(f">{sequence_id}\n{sequence}\n")
