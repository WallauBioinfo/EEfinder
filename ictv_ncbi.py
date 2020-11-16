#!/usr/bin/python3
# -*- coding: utf-8 -*-
from Bio import Entrez
import argparse, csv, re, time, subprocess, shlex
import xml.etree.ElementTree as ET
import json
import pandas as pd

parser = argparse.ArgumentParser(description = 'This script a csv file and returns protein information by viral family.')
parser.add_argument("-in", "--input", help="CSV file with 3 columns", required=True)
parser.add_argument("-p","--threads",help="Threads for cd-hit-est analyze.", default = 1, type=str)
parser.add_argument("-m","--memory",help="Memory in MB for cd-hit-est analyze",default = 2000, type=str)
args = parser.parse_args()
input_file = args.input
var_threads = args.threads
var_memory = args.memory

Entrez.email = 'zimmer.filipe@gmail.com'
Entrez.api_key = '29a8e7fd2add3e853ce3d3c18f140399a809'

def get_ids(response) -> list:
    j = json.loads(response.read())
    return list(j['esearchresult']['idlist'])

headers = ['Family','Genus','Species']
df = pd.read_csv(input_file, sep=',',header=None,names=headers)
for species in df.Species:
    specie_name = re.sub(r' ','_',species).rstrip('\n')
    try:   
        with open(specie_name+'.prot.fasta','a') as output_specie: 
            prids = get_ids(Entrez.esearch(db="Protein", term=F"{species}[Scientific Name]", retmode="json"))
            for prid in prids:
                fasta = Entrez.efetch(db="Protein", id=prid, rettype="fasta", retmode="text").read()
                output_specie.write(fasta)
                time.sleep(0.15)
            output_specie.close()
        with open(specie_name+'.prot.fasta.cd-hit.log','w') as cd_hit_log:
            cd_hit_out = specie_name+'.prot.fasta.cd'
            cd_hit_est_cmd = 'cd-hit -i '+specie_name+'.prot.fasta -o '+cd_hit_out+' -c 0.8 -aL 0.3 -g 1 -M '+str(var_memory)+' -T '+str(var_threads)+' -d 200'
            cd_hit_est_cmd = shlex.split(cd_hit_est_cmd)
            cd_hit_est_cmd_process = subprocess.Popen(cd_hit_est_cmd,stdout = cd_hit_log)
            cd_hit_est_cmd_process.wait()
        subprocess.call([F"sed -i -e 's/\[.*\]/{species}/g' {specie_name}.prot.fasta.cd"], shell=True)
        print(f'DONE for {species}')
    except:
        print(f'error for {species}')