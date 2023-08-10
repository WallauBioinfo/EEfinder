from Bio.Blast.Applications import NcbimakeblastdbCommandline
import subprocess
import shlex


def makeblastdb(data: str, db_type: str) -> None:
    clinedb = NcbimakeblastdbCommandline(dbtype=db_type, input_file=data)
    stdout, stderr = clinedb()


def makediamonddb(data: str, threads: int) -> None:
    clinedb = f"diamond makedb --db {data} --in {data} --threads {int(threads)} --matrix BLOSUM45"
    clinedb = shlex.split(clinedb)
    cmd_clinedb = subprocess.Popen(
        clinedb, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    cmd_clinedb.wait()


class MakeDB:
    """
    Creates BLAST or DIAMOND databases.

    Keyword arguments:
    mode: select the mode for mounting the database
    data: the database_file, parsed with -db argument
    db_type: 'nucl' or 'prot' strings
    threads: number of threads to be used in the throught the workflow
    """

    def __init__(self, mode: str, data: str, db_type: str, threads: int) -> object:
        self.mode = mode
        self.data = data
        self.db_type = db_type
        self.threads = threads

        self.make_db()

    def make_db(self) -> None:
        if self.mode in ["blastx", "tblastn"]:
            makeblastdb(self.data, self.db_type)
        else:
            makediamonddb(self.data, self.threads)
