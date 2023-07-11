from Bio.Blast.Applications import NcbimakeblastdbCommandline
import subprocess
import shlex


def makeblastdb(data, db_type):
    """
    This function creates blast databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    db_type - 'nucl' or 'prot' strings
    """
    clinedb = NcbimakeblastdbCommandline(dbtype=db_type, input_file=data)
    stdout, stderr = clinedb()


def makediamonddb(data, threads):
    """
    This function creates diamond databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    in_db_type - 'nucl' or 'prot' strings
    """
    clinedb = (
        f"diamond makedb --db {data} --in {data} --threads {threads} --matrix BLOSUM45"
    )
    clinedb = shlex.split(clinedb)
    cmd_clinedb = subprocess.Popen(clinedb, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cmd_clinedb.wait()


class MakeDB:
    def __init__(self, mode, data, db_type, threads):
        self.mode = mode
        self.data = data
        self.db_type = db_type
        self.threads = threads

        self.make_db()

    def make_db(self):
        if self.mode in ["blastx", "tblastn"]:
            makeblastdb(self.data, self.db_type)
        else:
            makediamonddb(self.data, self.threads)
