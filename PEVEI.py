#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Build in Python 3.6.5+

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://github.com/dezordi
"""

import pandas as pd
import argparse, csv, os, subprocess, shlex, sys, time, re
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from urllib.error import HTTPError
from socket import error as SocketError

#Setting arguments
parser = argparse.ArgumentParser(description = 'This script receives a fasta file, and predict regions of Endogenous Viral Elements')
parser.add_argument("-in", "--input", help="Fasta file", required=True)
parser.add_argument("-db", "--database", help="Viral proteins database file", required=True)
parser.add_argument("-ln", "--length", help="Length used to merge regions on bedtools merge, default = 100",type=int, default=100)
parser.add_argument("-p", "--threads", help="Threads for multi-thread analysis, default = 1",type=int, default=1)
parser.add_argument("-rm","--removetmp",help="Remove temporary files generated through analysis?, default = True", default = True, choices=['True','False'])
parser.add_argument("-mk","--makeblastdb",help="Make blast database?, default = True", default = True, choices=['True','False'])
#Storing argument on variables
args = parser.parse_args()
input_file = args.input
database_file = args.database
length_merge = args.length
threads = args.threads
temp_remove = args.removetmp
make_db = args.makeblastdb
#Connecting ncbi-email to Entrez package
#Entrez.email = input('Digite your e-mail to logging on NCBI: ') #add this line after...
#Functions:
def makeblastdb(data):
    """
    This function creates blast databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    """

    clinedb = NcbimakeblastdbCommandline(dbtype='prot', input_file = data)
    stdout, stderr = clinedb()
    return(print(f'Makeblastdb for {data}: DONE'))
def runblastx(query_file,database_file):
    """
    This function runs blastn.

    Keyword arguments:
    query_file - input_file - parsed with -in argument
    database - pre-formated blast database - parsed with -db argument
    """
    
    cline = NcbiblastxCommandline(query = query_file, db = database_file, out = query_file+'.blastx',
    outfmt = 6, evalue = 0.000001, num_threads = threads, matrix = 'BLOSUM45', max_intron_length = 90, soft_masking = 'true')
    stdout, stderr = cline()
    return(print(f'BLASTx with {query_file} against {database}: DONE'))

if __name__ == '__main__':
    """
    Main Routine
    This block of code is executed, whenever the script
    is started from the command line.
    """
    '''
    new_database = [] # list to recevei new records of database_file
    output_database = open(database_file+'.fmt','w')
    database_sequences = SeqIO.parse(database_file,"fasta")
    for record in database_sequences:
        record.id = re.sub(r'(?<=\.[0-9]).*$', '', record.id)
        new_database.append(record)
    #SeqIO.write(new_database,output_database,"fasta-2line")
    if make_db == True:
        makeblastdb(database_file+'.fmt')
    '''
    runblastx(input_file,database_file+'.fmt')    