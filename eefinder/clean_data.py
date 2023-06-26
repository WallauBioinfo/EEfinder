from Bio import SeqIO


def cut_seq(seq_file, cutoff, log):
    """
    This function remove contigs bellow the threshold.

    Keyword arguments:
    seq_file - fasta file, parsed with -in argument
    cutoff - cutoff length, parsed with -ln argument
    """
    new_sequences = []
    input_handle = open(seq_file, "r")
    output_handle = open(seq_file + ".fmt", "w")
    for record in SeqIO.parse(input_handle, "fasta"):
        if len(record.seq) >= int(cutoff):
            new_sequences.append(record)
    SeqIO.write(new_sequences, output_handle, "fasta")
    print(
        f"DONE: Sequences bellow than {cutoff} bp are removed from {seq_file} Results are stored in\n{seq_file}.fmt!",
        file=log,
    )
    return print(
        f"DONE: Sequences bellow than {cutoff} bp are removed from {seq_file} Results are stored in\n{seq_file}.fmt!"
    )


def mask_clean(in_file, m_per, log):
    """
    This function execute the cleanning step.

    Keyword arguments:
    in_file - fasta file, with putative EVEs.
    m_per - treshold masked percentage value, parsed with -mp argument.
    """

    sequences = {}
    for seq_record in SeqIO.parse(in_file, "fasta"):
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
        ) * 100 <= float(m_per):
            if sequence_id not in sequences:
                sequences[sequence_id] = sequence
    with open(in_file + ".cl", "w+") as output_file:
        for sequence_id, sequence in sequences.items():
            output_file.write(f">{sequence_id}\n{sequence}\n")
    print(
        f"DONE: Sequences with {m_per} percent of lower-case letters are removed, results are stored in {in_file}.cl!",
        file=log,
    )
    return print(
        f"DONE: Sequences with {m_per} percent of lower-case letters are removed, results are stored in {in_file}.cl!"
    )


class RemoveShort:
    def __init__(self, seq_file, cutoff, log):
        self.seq_file = seq_file
        self.cutoff = cutoff
        self.log = log

    def run_remove_short(self):
        cut_seq(self.seq_file, self.cutoff, self.log)


class MaskClean:
    def __init__(self, in_file, m_per, log):
        self.in_file = in_file
        self.m_per = m_per
        self.log = log

    def run_mask_clean(self):
        mask_clean(self.in_file, self.m_per, self.log)
