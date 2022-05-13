from Bio.Blast.Applications import NcbimakeblastdbCommandline
import subprocess
import shlex


def makeblastdb(data, db_type,make_db,log):
    """
    This function creates blast databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    db_type - 'nucl' or 'prot' strings
    """
    if make_db == True or db_type == 'nucl':
        clinedb = NcbimakeblastdbCommandline(dbtype=db_type, input_file=data)
        stdout, stderr = clinedb()
        print(f'DONE: Makeblastdb for {data}', file = log)
        return(print(f'DONE: Makeblastdb for {data}'))
    else:
        print(f'Pass Makeblastdb step', file = log)
        return(print(f'Pass Makeblastdb step.'))


def makediamonddb(data, threads,make_db,log):
    """
    This function creates diamond databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    in_db_type - 'nucl' or 'prot' strings
    """
    if make_db == True:
        clinedb = f"diamond makedb --db {data} --in {data} --threads {threads} --matrix BLOSUM45"
        clinedb = shlex.split(clinedb)
        cmd_clinedb = subprocess.Popen(clinedb)
        cmd_clinedb.wait()
        print(f'DONE: Makediamonddb for {data}', file = log)
        return(print(f'DONE: Makediamonddb for {data}'))
    else:
        print(f'Pass Makediamonddb step.', file = log)
        return(print(f'Pass Makediamonddb step.'))


class MakeDB():
    def __init__(self, mode, data, db_type, threads, make_db, log):
        self.mode = mode
        self.data = data
        self.db_type = db_type
        self.threads = threads
        self.make_db = make_db
        self.log = log

    def run_make_db(self):
        if self.mode in ['blastx', 'tblastn']:
            makeblastdb(self.data,self.db_type, self.make_db, self.log)
        else:
            makediamonddb(self.data, self.threads, self.make_db, self.log)