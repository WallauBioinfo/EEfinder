from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
import subprocess
import shlex


def runblastx(query_file: str, database_file: str, threads: int) -> None:
    cline = NcbiblastxCommandline(
        query=query_file,
        db=database_file,
        out=f"{query_file}.blastx",
        outfmt=6,
        word_size=3,
        evalue=0.00001,
        num_threads=threads,
        matrix="BLOSUM45",
        max_intron_length=100,
        soft_masking="true",
    )
    stdout, stderr = cline()


def rundiamond(query_file: str, database_file: str, threads: int, mode: str):
    clinedmd = f"diamond blastx -p {int(threads)} -d {database_file}.dmnd -f 6 -q {query_file} -o {query_file}.blastx -e 0.00001 --matrix BLOSUM45 -k 500 --max-hsps 0 --{mode}"
    clinedmd = shlex.split(clinedmd)
    cmd_clinedmd = subprocess.Popen(
        clinedmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    cmd_clinedmd.wait()


class SimilaritySearch:
    """
    Runs an BLAST or DIAMOND analysis

    Keyword arguments:

    query_file: file going to be used as query
    database_file: database created by make_db function
    threads: number of threads to be used in the throught the workflow
    mode: select the mode for run the similarity analisys
    """

    def __init__(
        self, query_file: str, database_file: str, threads: int, mode: str) -> object:
        self.query_file = query_file
        self.database_file = database_file
        self.threads = threads
        self.mode = mode

        self.similarity_search()

    def similarity_search(self) -> None:
        if self.mode == "blastx":
            runblastx(self.query_file, self.database_file, self.threads)
        else:
            rundiamond(self.query_file, self.database_file, self.threads, self.mode)
