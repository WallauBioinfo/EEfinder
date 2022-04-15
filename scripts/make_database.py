from Bio.Blast.Applications import NcbimakeblastdbCommandline
import subprocess, shlex

def makeblastdb(data,db_type):
    """
    This function creates blast databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    db_type - 'nucl' or 'prot' strings
    """
    clinedb = NcbimakeblastdbCommandline(dbtype=db_type, input_file = data)
    stdout, stderr = clinedb()
    return(print(f'DONE: Makeblastdb for {data}\n'))

def makediamonddb(data, threads):
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
    return(print(f'DONE: Makediamonddb for {data}\n'))

class MakeDB():
    def __init__(self,mode,data,db_type,threads):
        self.mode = mode
        self.data = data
        self.db_type = db_type
        self.threads = threads
    
    def run_make_db(self):
        if self.mode == 'blastx':
            makeblastdb(data = self.data, db_type = self.db_type)
            #print(f"makeblastdb(data = {self.data}, db_type = {self.db_type})")
        else:
            makediamonddb(data = self.data, threads = self.threads)
            #print(f"makediamonddb(data = {self.data}, threads = {self.threads})")