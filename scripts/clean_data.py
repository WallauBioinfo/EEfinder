from Bio import SeqIO

def cut_seq(seq_file,cutoff):
    """
    This function remove contigs bellow the threshold.

    Keyword arguments:
    seq_file - fasta file, parsed with -in argument
    cutoff - cutoff length, parsed with -ln argument
    """
    new_sequences = []
    input_handle=open(seq_file,'r')
    output_handle=open(seq_file+'.fmt','w')
    for record in SeqIO.parse(input_handle,'fasta'):
        if len(record.seq) >= int(cutoff):
            new_sequences.append(record)
    SeqIO.write(new_sequences,output_handle,"fasta")
    return(print(f'DONE: Sequences bellow than {cutoff} bp are removed from {seq_file},\nresults are stored in {seq_file}.fmt!\n'))

def mask_clean(in_file, m_per):
    """
    This function execute the cleanning step.

    Keyword arguments:
    in_file - fasta file, with putative EVEs.
    m_per - treshold masked percentage value, parsed with -mp argument.
    """
    sequences={}
    for seq_record in SeqIO.parse(in_file, "fasta"):
        sequence = str(seq_record.seq)
        if (float(sequence.count("a")+sequence.count("t")+sequence.count("c")+sequence.count("g")+sequence.count("n")+sequence.count("N"))/float(len(sequence)))*100 <= float(m_per):
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
    with open(in_file+'.cl', "w+") as output_file:
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
    return(print(f'Sequences with {m_per} percent of lower-case letters are removed, results are stored in {in_file}.clean!'))

class RemoveShort():
    def __init__(self, seq_file, cutoff):
        self.seq_file = seq_file
        self.cutoff = cutoff
    
    def run_remove_short(self):
        cut_seq(self.seq_file, self.cutoff)

class MaskClean():
    def __init__(self, in_file, m_per):
        self.in_file = in_file
        self.cutoff = cutoff

    def run_mask_clean(self):
        mask_clean(self.in_file, self.m_per)
