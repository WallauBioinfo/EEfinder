from Bio.Blast.Applications import NcbimakeblastdbCommandline
import subprocess
import shlex


def makeblastdb(data, db_type, log):
    """
    This function creates blast databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    db_type - 'nucl' or 'prot' strings
    """
    clinedb = NcbimakeblastdbCommandline(dbtype=db_type, input_file=data)
    stdout, stderr = clinedb()
    print(f"DONE: Makeblastdb for {data}", file=log)
    return print(f"DONE: Makeblastdb for {data}")

def makediamonddb(data, threads, log):
    """
    This function creates diamond databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    in_db_type - 'nucl' or 'prot' strings
    """
    clinedb = f"diamond makedb --db {data} --in {data} --threads {threads} --matrix BLOSUM45"
    clinedb = shlex.split(clinedb)
    cmd_clinedb = subprocess.Popen(clinedb)
    cmd_clinedb.wait()
    print(f"DONE: Makediamonddb for {data}", file=log)
    return print(f"DONE: Makediamonddb for {data}")


class MakeDB:
    def __init__(self, mode, data, db_type, threads, log):
        self.mode = mode
        self.data = data
        self.db_type = db_type
        self.threads = threads
        self.log = log

    def run_make_db(self):
        if self.mode in ["blastx", "tblastn"]:
            makeblastdb(self.data, self.db_type, self.log)
        else:
            makediamonddb(self.data, self.threads, self.log)
