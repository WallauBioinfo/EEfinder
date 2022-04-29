from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
import subprocess
import shlex


def runblastx(query_file, database_file, threads, log):
    """
    This function runs Blastx analisis.

    Keyword arguments:
    query_file - input_file - parsed with -in argument
    database - pre-formated blast database - parsed with -db argument

    Blast arguments:
    max_intron_length = "Length of the largest intron allowed in a translated nucleotide sequence when linking multiple distinct alignments (a negative value disables linking)."
    soft_masking = "Apply filtering locations as soft masks (i.e., only for finding initial matches)."
    """

    cline = NcbiblastxCommandline(query=query_file, db=database_file, out=f"{query_file}.blastx",
                                  outfmt=6, word_size=3, evalue=0.00001, num_threads=threads, matrix='BLOSUM45', max_intron_length=100, soft_masking='true')
    stdout, stderr = cline()
    print(f'DONE: BLASTx with {query_file}\nagainst {database_file}',file = log)
    return(print(f'DONE: BLASTx with {query_file}\nagainst {database_file}'))


def runtblastn(query_file, database_file, threads, output, log):  # colocar opções de blastn ou blastx
    """
    This function runs blast with putative EVEs against host genes database.

    Keyword arguments:
    query_file - input_file - parsed with -in argument
    database - pre-formated blast database - parsed with -db argument
    """

    cline = NcbitblastnCommandline(query=query_file, db=database_file, out=output,
                                   outfmt=6, evalue=0.00001, num_threads=threads,  max_intron_length=100, soft_masking='true')
    stdout, stderr = cline()
    print(f'DONE: tBLASTn with {query_file} against {database_file}', file = log)
    return(print(f'DONE: tBLASTn with {query_file} against {database_file}'))


def rundiamond(query_file, database_file, threads, mode, log):
    """
    This function runs Diamond analysis.

    Keyword arguments:
    query_file - input_file - parsed with -in argument
    database - pre-formated database - parsed with -db argument

    Diamond arguments:
    -k 500 --max-hsps 0 are the correspondents to  max_target_seqs 500 and max_hsps None in BLAST parameters.
    the BLAST arguments word_size, max_intron_length and soft_masking are not available in DIAMOND tool.
    """

    clinedmd = f"diamond blastx -p {threads} -d {database_file}.dmnd -f 6 -q {query_file} -o {query_file}.blastx -e 0.00001 --matrix BLOSUM45 -k 500 --max-hsps 0 --{mode}"
    clinedmd = shlex.split(clinedmd)
    cmd_clinedmd = subprocess.Popen(clinedmd)
    cmd_clinedmd.wait()
    print(f'DONE: DIAMOND {mode} analysis with {query_file} against {database_file}', file = log)
    return(print(f'DONE: DIAMOND {mode} analysis with {query_file} against {database_file}'))


class SimilaritySearch():
    def __init__(self, query_file, database_file, threads, mode, log):
        self.query_file = query_file
        self.database_file = database_file
        self.threads = threads
        self.mode = mode
        self.log = log

    def run_similarity_search(self):
        if self.mode == 'blastx':
            runblastx(self.query_file, self.database_file, self.threads,self.log)
        else:
            rundiamond(self.query_file, self.database_file,
                       self.threads, self.mode,self.log)


class FlankSearch():
    def __init__(self, query_file, database_file, threads, output, log):
        self.query_file = query_file
        self.database_file = database_file
        self.threads = threads
        self.output = output
        self.log = log

    def run_flank_search(self):
        runtblastn(self.query_file, self.database_file,
                   self.threads, self.output, self.log)
